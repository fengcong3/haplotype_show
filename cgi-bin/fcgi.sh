#./fcgi.py --method=prefork/threaded minspare=50 maxspare=50 maxchildren=1000 
spawn-fcgi  -f /var/www/g3rp/cgi-bin/fcgi.py -a 127.0.0.1 -p 8008  -F 5 -P /var/www/g3rp/cgi-bin/fcgi.pid -u ghost1
