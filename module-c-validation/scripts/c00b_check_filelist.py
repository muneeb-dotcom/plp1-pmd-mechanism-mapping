import requests
r = requests.get("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277705/suppl/filelist.txt", timeout=60)
print(r.text)