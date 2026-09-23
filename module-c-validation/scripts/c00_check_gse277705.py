import requests

url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277705/suppl/"
r = requests.get(url, timeout=60)
r.raise_for_status()
print(r.text)