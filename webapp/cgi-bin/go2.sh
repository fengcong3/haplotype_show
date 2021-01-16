#./fcgi.py --method=prefork/threaded minspare=50 maxspare=50 maxchildren=1000 
spawn-fcgi  -f /var/www/hap.bioinf.club/webapp/cgi-bin/down.py -a 127.0.0.1 -p 8011  -F 5 -P /var/www/hap.bioinf.club/webapp/cgi-bin/fcgi2.pid -u ghost1
