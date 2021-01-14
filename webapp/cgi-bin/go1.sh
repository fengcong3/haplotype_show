#./fcgi.py --method=prefork/threaded minspare=50 maxspare=50 maxchildren=1000 
spawn-fcgi  -f /var/www/hap.bioinf.club/webapp/cgi-bin/toiTol.py -a 127.0.0.1 -p 8010  -F 5 -P /var/www/hap.bioinf.club/webapp/cgi-bin/fcgi1.pid -u ghost1
