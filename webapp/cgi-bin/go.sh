#./fcgi.py --method=prefork/threaded minspare=50 maxspare=50 maxchildren=1000 
spawn-fcgi  -f /var/www/hap.bioinf.club/webapp/cgi-bin/hap-show.py -a 127.0.0.1 -p 8009  -F 5 -P /var/www/hap.bioinf.club/webapp/cgi-bin/fcgi.pid -u ghost1
