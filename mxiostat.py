#!/usr/bin/python
#
# A version of iostat that actually reports real numbers taken straight
# from the kernel. This uses the 2.6 /proc/diskstats stuff.
#
# Usage: mxiostat.py [-qvaH] [-D DEV,...] [-c COUNT] [[DEV|FIELD,...] [DELAY]]
#
# -q: don't print the periodic header.
# -v: be more verbose in some situations.
# -c COUNT: exit after printing COUNT reports.
# -a: report on particular stat(s) for all disks, instead of all stats for
#	one disk
# -C: compact reporting for multiple fields in -a mode; all fields for all
#     disks are printed on one line.
# -V: 'vertical' report for multiple disks in -a mode, where each disk is
#     printed on one line
# -D DEV,...: in -a, report on only these device(s) instead of all devices
# -X DEV,...: in -a, exclude these devices from the report
#
# DEV is the device to report stats on; it defaults to 'sda'.
# FIELD,... is the field or fields to report on for all disks; it defaults
#	to 'util'. 'all' may be used for all fields, and it is the default
#	for -V.
# DELAY is the amount of time to delay between report lines, in
# seconds; it defaults to one (1) second. You can use fractional
# seconds if you want to.
#
# -a has three display modes: fields down the side and disks horizontally
# (the default), disks down the side and fields horizontally (-V), and a
# compact form where everything is printed on one line (-C).
#
# All statistics are per-second ones, regardless of the delay, except for
# 'act' which is an instantaneous number (how many IOs are in flight RIGHT
# NOW, as the program reads the data).
# 'rwait', 'wwait', and 'await' are in milliseconds. 'util' is a percentage.
# 'rsect' and 'wsect' are 512-byte units. (We could print rkB and wkB, but
# the author wants this to be as close to raw as possible.)
#
# See the end of the code for more information on fields.

# This avoids all sorts of irritating 'float(...)/...' usage and lets
# us just divide straight.
from __future__ import division

import sys, time, getopt, os.path
import re

# This is not a user-servicable part. It controls whether we barf up
# messages when delta-computed fields go negative.
debug = 0

class IOEx(Exception):
	pass

# Resolve a filesystem mount point into a device name suitable for
# finding in /proc/partitions et al.
# Note that in /proc/mounts, loopback mounts show as being mounted
# on the real device. It's only in 'mount' that they come up funny.
def find_device(fsname):
	fp = open("/proc/mounts", "r")
	for ln in fp:
		flds = ln.strip().split()
		if flds[1] != fsname:
			continue
		if flds[2] not in ('ext2', 'ext3', 'vfat', ):
			raise IOEx, "cannot deal with filesystem type %s for %s" % (flds[2], fsname)
		return flds[0]
	raise IOEx, "could not find /proc/mounts entry for filesystem %s" % fsname
lvre = re.compile("\[[A-Z ]+\] ([^\s]+)\s+\d+ /\d+ ")
def resolve_lvm(devname):
	# Split LVM name into VG and LV.
	ab = devname.split("/")
	fp = open("/proc/lvm/global", "r")
	curvg = ""; physdev = ""
	inlv = 0
	for ln in fp:
		ln = ln.strip()
		if ln.startswith("VG: "):
			curvg = ln.split()[1]
			inlv = 0
		elif ln.startswith("PV: "):
			physdev = ln.split()[2]
		elif ln.startswith("LVs: ") or ln.startswith("LV: "):
			inlv = 1
			ln = ln.split(None, 1)[1]
		if inlv:
			mo = lvre.match(ln)
			if not mo:
				inlv = 0
				continue
			if mo.group(1) == ab[1] and curvg == ab[0]:
				return physdev
	raise IOEx, "Cannot find the physical volume for LVM volume /dev/%s" % devname	
def resolve_fsname(fsname):
	dn = find_device(fsname)
	if dn.startswith("/dev/md"):
		raise IOEx, "filesystem %s is on MD device %s, we cannot deal with those" % (fsname, dn)
	elif not dn.startswith("/dev/"):
		raise IOEx, "filesystem %s is on %s, which I cannot deal with" % (fsname, dn)
	# zap off the /dev.
	dn = dn[len("/dev/"):]
	# Okay, are we on an LVM volume? (This is an ugly heuristic.)
	if not '/' in dn:
		return dn
	return resolve_lvm(dn)

# Subtract two unsigned 32-bit integers in the presence of rollover.
# The best we can do is assume that the rollover only happens once.
#
# While the theory that negative number freakouts are caused by
# (unsigned) integer rollovers is a nice one and attractive and all
# that, it turns out to be completely wrong. While they are possible,
# dumping the numbers involved makes it clear that they are not the
# actual cause; the 'b' number is only around 2**25 to 2**26 or so.
# At the moment, the author is inclined to blame (lack of) locking on
# the gendisk stats, since the stat update caused by reading them from
# /proc/partitions can likely race with stat updates from interrupt
# level.
#
# Still, we preserve this code for posterity just in case.
maxU = 2**32-1
def uintSub(a, b):
	# First: some of these numbers are not actually merely
	# unsigned integers; some of them are long longs. Cope
	# with this case. (I *think* this doesn't happen, but
	# better safe than sorry.)
	if a > maxU or b > maxU or a >= b:
		return a - b
	res = (a - b) % maxU
	#print "uintSub:", a, b, res
	return res

# This stores the IO statistics (for a single device object) harvested
# from /proc/partitions. We define subtraction on iostat objects as a
# delta operation, which means we leave 'running' intact.
class iostats:
	fields = ('rio', 'rmerge', 'rsect', 'ruse', 'wio', 'wmerge', 'wsect',
		  'wuse', 'running', 'use', 'aveq', 'rduse', 'wduse',
		  'rtime', 'wtime', 'raveq', 'waveq', )
	def __init__(self, lst = None):
		if not lst:
			for f in self.fields:
				setattr(self, f, 0)
		else:
			for x in xrange(0, min(len(self.fields), len(lst))):
				setattr(self, self.fields[x], lst[x])
			# None is the special 'not supported' value.
			for x in xrange(len(lst), len(self.fields)):
				setattr(self, self.fields[x], None)
	def __sub__(self, other):
		n = iostats()
		for f in self.fields:
			# 'running' is an instantaneous number, so we just
			# copy it, not subtract it.
			if f == 'running':
				setattr(n, f, getattr(self, f))
			# propagate 'not supported' magically.
			elif getattr(self, f) is None:
				setattr(n, f, None)
			else:
				#res = uintSub(getattr(self, f),
				#	      getattr(other, f))
				res = getattr(self, f) - getattr(other, f)
				if res < 0:
					if debug:
						sys.stderr.write("xiostat: %s freakout: %d - %d = %d\n" % (f, getattr(self, f), getattr(other, f), res))
					res = -1
				setattr(n, f, res)
		return n
	def __str__(self):
		return "<iostat: %s>" % ", ".join(["%s: %d" % (f, getattr(self, f)) for f in self.fields])

# Read /proc/uptimes.
# We return the uptime in seconds (the first field). Note that this can
# have decimal digits, so we must return it as a floating point number.
def getuptime(fn):
	fp = open(fn, "r")
	l = fp.readline()
	fp.close()
	return float(l.split()[0])

# Read /proc/partitions, returning an iostat object.
# This will alternately read /proc/diskstats, returning the same.
# We tell the difference based on black magic.
def getdiskstats(fn):
	fp = open(fn, "r")
	devAt = None; lineLen = None
	rdevs = {}
	for l in fp:
		n = l.strip().split()
		# first line: determine format.
		if devAt is None:
			if n[0] == "major":
				devAt = 3; lineLen = 15
				continue
			else:
				devAt = 2; lineLen = 14
		# ugly hack to screen out everything except real
		# disks.
		devname = n[devAt]
		if len(n) < lineLen:
			continue
		rdevs[devname] = iostats([long(x) for x in n[devAt+1:]])
	fp.close()
	return rdevs
def getdiskstat(fn, dev):
	d = getdiskstats(fn)
	return d.get(dev, None)

# Out of the universe of things with valid-looking disk stats, produce
# a list of things that are real disks.
def getrealdisks(dkstats):
	dlist = []
	for devn in dkstats.keys():
		# gets partitions, loopN, md devices, etc etc
		if devn[-1].isdigit():
			continue
		if devn == 'hda':
			continue
		dlist.append(devn)
	dlist.sort()
	return dlist

# Return a decimalized string that will fit into a field that is fw wide.
# We never use more than two decimal places, but we can use one or
# zero depending on the size of 'num'.
# Life is complicated by our need to cope with the rounding that goes
# on with floating point numbers during printing.
def decimalize(num, fw):
	maxtwo = 10 ** (fw-3)
	maxone = 10 ** (fw-2)
	# Sometimes formatting a floating point number for printing with
	# a limited number of digits can make it larger, causing us to
	# potentially overflow a field width. We use round() beforehand
	# here to cope with that; otherwise we would get periodic results
	# of '100.0' for a four-digit field, for example, when we get
	# handed '99.999' as the number to print.
	if round(num, 1) >= maxone or long(num) == num:
		# '%d' apparently truncates (?!); oh well, cope.
		return ("%%%dd" % fw) % round(num, 0)
	elif round(num, 2) >= maxtwo:
		return ("%%%d.1f" % fw) % num
	else:
		return ("%%%d.2f" % fw) % num

# This calculates the numeric value from the field name, the IO delta
# structure, and the time delta. (The time delta is in seconds.)  Some
# fields are more or less straight from the iod (sometimes recomputed
# to per-second values); others are calculated from multiple others.
persecs = ('rio', 'rmerge', 'rsect', 'wio', 'wmerge', 'wsect',)
def calcval(iod, td, field):
	# if b is nonzero, return a/b; otherwise return 0.
	def zerodiv(a, b):
		# In rare cases (iod.use especially) the denominator can
		# be negative instead of the numerator.
		if a is None:	return None
		elif a < 0:	return a
		elif b < 0:	return b
		if b:	return a / b
		else:	return 0
	if field == 'act':		return iod.running
	elif field == 'rwait':		return zerodiv(iod.ruse, iod.rio)
	elif field == 'wwait':		return zerodiv(iod.wuse, iod.wio)
	elif field == 'await':
		# await is invalid if either of the components is invalid,
		# so we must ripple this invalidity through.
		if iod.ruse == -1 or iod.wuse == -1:
			return -1
		return zerodiv(iod.ruse + iod.wuse, iod.rio + iod.wio)
	elif field == 'rdwait':		return zerodiv(iod.rduse, iod.rio)
	elif field == 'wdwait':		return zerodiv(iod.wduse, iod.wio)
	elif field == 'adwait':
		if iod.rduse is None or iod.wduse is None:
			return None
		if iod.rduse == -1 or iod.wduse == -1:
			return -1
		return zerodiv(iod.rduse + iod.wduse, iod.rio + iod.wio)
	elif field == 'aveq':		return zerodiv(iod.aveq, iod.use)
	elif field == 'raveq':		return zerodiv(iod.raveq, iod.rtime)
	elif field == 'waveq':		return zerodiv(iod.waveq, iod.wtime)
	elif field == 'rgrp':		return zerodiv(iod.rsect, iod.rio)
	elif field == 'wgrp':		return zerodiv(iod.wsect, iod.wio)
	elif field == 'agrp':
		return zerodiv(iod.rsect + iod.wsect, iod.rio + iod.wio)
	# rolled up sums
	elif field == 'aio':		return (iod.rio + iod.wio) / td
	elif field == 'asect':		return (iod.rsect + iod.wsect) / td
	# 'util' is a percentage, and 'use' is in milliseconds.
	# This makes it work out nicely.
	elif field == 'util':
		# 'use' can apparently go negative too. Who knew? Sigh.
		if iod.use < 0:
			return -1
		return iod.use / (td*10)
	elif field == 'rutil':
		if iod.rtime is None:
			return None
		if iod.rtime < 0:
			return -1
		return iod.rtime / (td*10)
	elif field == 'wutil':
		if iod.wtime is None:
			return None
		if iod.wtime < 0:
			return -1
		return iod.wtime / (td*10)
	elif hasattr(iod, field) and field in persecs:
		return getattr(iod, field) / td
	elif hasattr(iod, field):
		return getattr(iod, field)
	else:
		raise IOEx, "don't know how to compute %s" % field

# The field width of each field.
fwidth = {'act': 4, 'rio': 6, 'wio': 6, 'rmerge': 6, 'wmerge': 6,
	  'rsect': 6, 'wsect': 6, 'rwait': 5, 'wwait': 5, 'await': 5,
	  'rgrp': 5, 'wgrp': 5, 'agrp': 5,
	  'rdwait': 6, 'wdwait': 6, 'adwait': 6,
	  'raveq': 5, 'waveq': 5, 'rutil': 5, 'wutil': 5,
	  'aveq': 4, 'util': 4,
	  'aio': 6, 'asect': 6, }
# The order in which fields are displayed in the output.
forder = ('act',
	  'rio', 'rmerge', 'rsect', 'rwait', 'rgrp',
	  'wio', 'wmerge', 'wsect', 'wwait', 'wgrp',
	  'agrp', 'aveq', 'await', 'util',
	  'rdwait', 'wdwait', 'adwait', 'rutil', 'wutil', 'raveq', 'waveq', )

# Display an iostat delta, using the given time delta to convert the
# relevant numbers to a per-second count. The time delta is in seconds.
def display(iod, td):
	outl = []
	for fn in forder:
		r = calcval(iod, td, fn)
		if r is None:
			continue
		outl.append(decimalize(r, fwidth[fn]))
	print " ".join(outl)
	# And flush our output if we are writing to a file or a pipe for
	# logging purposes:
	sys.stdout.flush()
# Print the report header.
def header(iod, td):
	outl = []
	for fn in forder:
		r = calcval(iod, td, fn)
		if r is None:
			continue
		outl.append(("%%%ds" % fwidth[fn]) % fn)
	print " ".join(outl)

# Print the mass report header
def fmtto(txt, width):
	return ("%%%ds" % width) % txt
def mheader(fields, devl):
	m1 = max([len(x) for x in devl])
	m2 = max([max(len(x), fwidth[x]) for x in fields])
	ml = max(m1, m2)
	outl = []
	if len(fields) == 1:
		flist = [fields[0],] + devl
	else:
		flist = ["",] + devl
	for fn in flist:
		outl.append(fmtto(fn, ml))
	print " ".join(outl)
	return ml

# Print the mass report line.
def mdisplay(field, deltas, td, mw, fname):
	outl = [fmtto(fname, mw)]
	for d in deltas:
		r = calcval(d, td, field)
		outl.append(decimalize(r, mw))
	print " ".join(outl)
	sys.stdout.flush()

# 'horizontal' display of multiple fields
# field widths vary across things.
def hwidth(fieldnames, devl):
	devwidth = max([len(x) for x in devl])
	fres = []
	for field in fieldnames:
		mw = max(devwidth, fwidth[field], len(field))
		fres.append((field, mw))
	return fres
def hheader(fields, devl):
	outl = []
	for fn, fw in fields:
		outl.append(fmtto(fn, fw))
		for d in devl:
			outl.append(fmtto(d, fw))
	print " ".join(outl)
def hdisplay(fields, deltas, td):
	outl = []
	for fn, fw in fields:
		outl.append(fmtto("", fw))
		for d in deltas:
			r = calcval(d, td, fn)
			outl.append(decimalize(r, fw))
	print " ".join(outl)
	sys.stdout.flush()

# 'vertical' display of multiple values.
# NNGH HEAD HURTS HATE IT
def vdisplay(dt, dname, td, flist):
	outl = []
	outl.append("%-6s" % dname)
	for fn in flist:
		r = calcval(dt, td, fn)
		if r is None:
			continue
		outl.append(decimalize(r, fwidth[fn]))
	print " ".join(outl)
	# And flush our output if we are writing to a file or a pipe for
	# logging purposes:
	sys.stdout.flush()
# Print the report header.
def vheader(dt, flist):
	outl = ["%-6s" % " ",]
	for fn in flist:
		r = calcval(dt, 1, fn)
		if r is None:
			continue
		outl.append(("%%%ds" % fwidth[fn]) % fn)
	print " ".join(outl)

# Our main processing loop proceeds by getting the initial stats,
# then looping around sleeping the delay time, getting new stats,
# computing the delta, and finally displaying them.
def statloop(every, dev, showheader = 1, max = 0):
	oldUt = getuptime("/proc/uptime")
	if os.path.exists("/proc/diskstats"):
		statFile = "/proc/diskstats"
	else:
		statFile = "/proc/partitions"
	oldSt = getdiskstat(statFile, dev)
	if not oldSt:
		raise IOEx, "cannot get starting stats for %s" % dev
	# We produce the header immediately as a reassurance that we
	# are actually doing something, since it could be some time
	# before the first display (since it comes after an 'every'
	# interval).
	if showheader:
		header(oldSt, 1)
	lc = 1
	while (max == 0) or lc < max:
		time.sleep(every)
		newUt = getuptime("/proc/uptime")
		newSt = getdiskstat(statFile, dev)
		if not newSt:
			raise IOEx, "cannot get new stats for %s" % dev
		td = newUt - oldUt
		iosd = newSt - oldSt
		# If there is no time delta, we punt.
		# We compute the time delta explicitly (from /proc/uptime's
		# uptime time) because on a loaded system we may have slept
		# for (much) longer than 'every' seconds.
		if td == 0:
			continue
		if showheader and lc % 22 == 0:
			header(iosd, td)
		display(iosd, td)
		oldSt = newSt
		oldUt = newUt
		lc += 1

# Our main processing loop proceeds by getting the initial stats,
# then looping around sleeping the delay time, getting new stats,
# computing the delta, and finally displaying them.
def massloop(every, fields, devlist, showheader, max, disp, excl):
	oldUt = getuptime("/proc/uptime")
	if os.path.exists("/proc/diskstats"):
		statFile = "/proc/diskstats"
	else:
		statFile = "/proc/partitions"
	oldSts = getdiskstats(statFile)
	if not oldSts:
		raise IOEx, "cannot get starting stats"
	odevs = set(getrealdisks(oldSts))
	devl = list(odevs)
	devl.sort()
	if devlist:
		ns = list(set(devlist) - set(oldSts.keys()))
		ns.sort()
		if ns:
			raise IOEx, "cannot get starting stats for: " + \
			      " ".join(ns)
	else:
		devlist = devl
	if excl:
		for e in excl:
			if e in devlist:
				devlist.remove(e)

	# evil hack
	if disp == 'horiz':
		flist = hwidth(fields, devlist)
		
	# We produce the header immediately as a reassurance that we
	# are actually doing something, since it could be some time
	# before the first display (since it comes after an 'every'
	# interval).
	if showheader:
		if disp == 'horiz':
			hheader(flist, devlist)
		elif disp == 'vert':
			vheader(oldSts[devlist[0]], fields)
		else:
			mw = mheader(fields, devlist)
	lc = 1
	lines = 1
	while (max == 0) or lc < max:
		time.sleep(every)
		newUt = getuptime("/proc/uptime")
		newSts = getdiskstats(statFile)
		if not newSts:
			raise IOEx, "cannot get new stats"
		# Check that we have the same disks.
		ndevs = set(getrealdisks(newSts))
		if ndevs != odevs:
			raise IOEx, "disks appeared or disappeared"
		td = newUt - oldUt

		# If there is no time delta, we punt.
		# We compute the time delta explicitly (from /proc/uptime's
		# uptime time) because on a loaded system we may have slept
		# for (much) longer than 'every' seconds.
		if td == 0:
			continue

		# Compute the difference for every device, now that we
		# know that there is a time difference to display results
		# for.
		deltas = [newSts[x] - oldSts[x] for x in devlist]

		# Show stuff.
		if showheader and lines >= 22:
			if disp == 'horiz':
				hheader(flist, devlist)
			elif disp == 'vert':
				vheader(deltas[0], fields)
			else:
				mheader(fields, devlist)
			lines = 0
		if disp == 'horiz':
			hdisplay(flist, deltas, td)
			lines += 1
		elif disp == 'vert':
			for i in range(len(devlist)):
				vdisplay(deltas[i], devlist[i], td, fields)
				lines += 1
			# force a header line for every report if we have
			# enough devices that two reports always overflow.
			if len(devlist) >= 11:
				lines = 22
		else:
			for field in fields:
				if len(fields) == 1:
					fname = ""
				else:
					fname = field
				mdisplay(field, deltas, td, mw, fname)
				lines += 1
			# force a header line every report as above (but
			# for fields).
			if len(fields) >= 11:
				lines = 22
		oldSts = newSts
		oldUt = newUt
		lc += 1

def error(msg):
	sys.stderr.write("%s: %s\n" % (sys.argv[0], msg))
def die(msg):
	error(msg)
	sys.exit(1)
def usage():
	sys.stderr.write("usage: mxiostat [-avqCV] [-D disk,disk,...] [-X disk,disk,...] [-c COUNT] [[DEV|FIELD[,FIELD]] [DELAY]]\n")
	sys.stderr.write("\tDEV defaults to sda, FIELD to util, DELAY to 1 second\n")
	sys.exit(1)
def process(args):
	repmax = 0; every = 1; dev = "sda"; showheader = 1
	verbose = 0; alldevs = 0; horiz = 0; vert = 0
	devlist = []; excl = []
	try:
		opt, arg = getopt.getopt(args, 'aqc:vD:CVX:', [])
	except getopt.GetoptError, e:
		error(str(e))
		usage()
	for o, a in opt:
		if o == '-q':
			showheader = 0
		elif o == '-a':
			alldevs = 1
			fields = ['util',]
		elif o == '-C':
			horiz = 1
		elif o == '-V':
			vert = 1
			fields = forder
		elif o == '-v':
			verbose += 1
		elif o == '-D':
			devlist = a.split(',')
		elif o == '-X':
			excl = a.split(',')
		elif o == '-c':
			try:
				repmax = int(a)
			except ValueError:
				die("-c argument '%s' is not an integer" % a)
		else:
			die("Chris failed to handle switch '%s'" % o)
	if len(arg) > 2:
		error("Extra arguments supplied.")
		usage()
	if len(arg) == 2:
		try:
			every = float(arg[1])
		except ValueError:
			error("Cannot turn delay '%s' into a number" % arg[1])
			usage()
	if len(arg) > 0 and not alldevs:
		dev = arg[0]
		if dev[0] == "/":
			dev = resolve_fsname(dev)
			if verbose:
				print "Displaying statistics for", dev
	elif len(arg) > 0:
		# field name(s)
		if arg[0] == 'all':
			fields = forder
		else:
			fields = arg[0].split(',')
			for f in fields:
				if f not in fwidth:
					die("no field %s" % f)
	if not (horiz or vert):
		display = "normal"
	elif horiz and vert:
		die("-H and -V are incompatible")
	elif horiz:
		display = "horiz"
	elif vert:
		display = "vert"
	try:
		if alldevs:
			massloop(every, fields, devlist, showheader, \
				 repmax, display, excl)
		else:
			statloop(every, dev, showheader, repmax)
		sys.exit(0)
	except IOEx, e:
		error(str(e))
		sys.exit(1)
	except EnvironmentError, e:
		error("System error: %s" % str(e))
		sys.exit(1)

if __name__ == "__main__":
	process(sys.argv[1:])

#
# -----
# Field values and what they mean:
#
# rio		number of read IO requests completed
# rmerge	number of submitted read requests that were merged into
#		existing requests.
# rsect		number of read IO sectors submitted (always 512-byte sectors)
# rgrp		rsect / rio, average number of sectors per request
#		(more or less)
# rwait		average wait time for read requests *in milliseconds*
#
# w*		as with the r* versions, but for writes
#
# agrp		overall average number of sectors per request
# await		overall average wait time for requests (in msec)
# aio		overall requests submitted (rio + wio)
# asect		overall sector count (rsect + wsect)
#
# act		number of IO requests in the queue at this instant
# aveq		average number of requests in the queue
# util		percent utilization of the drive (ie, % busy)
#
# Most of these fields are whatever-per-second; act is an instantaneous
# number. See also http://utcc.utoronto.ca/~cks/space/blog/linux/DiskIOStats
# for some additional notes and exact details.
#
# (additional fields visible in the source are vestigial)
