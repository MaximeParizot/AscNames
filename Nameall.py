import os  

filelist=os.listdir('OrthologyRelationships')
d={}
d['1ToMany']=filelist[0]
d['1To1']=filelist[2]
d['manToMany']=filelist[3]

for file in d.keys(): 
    commande='python3 Nametunicates.py -f OrthologyRelationships/'+d[file]+' -out Outputs/'+file+' -lim 100 -check n'
    os.system(commande)