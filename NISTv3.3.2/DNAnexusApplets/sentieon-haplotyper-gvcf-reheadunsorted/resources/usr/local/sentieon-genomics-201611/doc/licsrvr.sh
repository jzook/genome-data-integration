#!/bin/sh
### BEGIN INIT INFO
# Provides: sentieon-licsrvr
# Required-Start: $local_fs $network
# Required-Stop:  $local_fs $network
# Default-Start:  2 3 4 5
# Default-Stop: 0 1 6
# Short-Description: start and stop Sentieon license server
# Description: Sentieon license server provides floating license support for
#	Sentieon bioinformatic tools.
### END INIT INFO

# source function library
. /lib/lsb/init-functions

# pull in sysconfig settings
[ -f /etc/sysconfig/licsrvr ] && . /etc/sysconfig/licsrvr
[ -f /etc/default/licsrvr ] && . /etc/default/licsrvr

RETVAL=0
prog="licsrvr"

licsrvr=${licsrvr:-/usr/sbin/licsrvr}
licfile=${licfile:-/etc/sentieon.lic}
logfile=${logfile:-/var/log/licsrvr.log}

start() {
	echo -n "Starting $prog:"
	$licsrvr --log $logfile --start $licfile
	RETVAL=$?
	echo
}

stop() {
	echo -n "Stopping $prog:"
	$licsrvr --stop $licfile
	RETVAL=$?
	echo
}

status() {
	pid=$(pidofproc $prog)
	RETVAL=$?
	case "$RETVAL" in
	0)	echo "$prog (pid $pid) is running...";;
	3)	echo "$prog is stopped";;
	4)	echo "$prog status unknown";;
	esac
}

case "$1" in
	start)
		start
		;;
	stop)
		stop
		;;
	restart)
		stop
		sleep 3
		start
		;;
	status)
		status
		;;
	*)
		echo "Usage: $0 {start|stop|restart|status}"
		RETVAL=1
esac
exit $RETVAL
