import os  

filelist=os.listdir('OrthologyRelationships')
d={}
d['1ToMany']='rootedTre_orthology1ToMany.txt'
d['1To1']='rootedTre_orthology1To1.txt'
d['manyToMany']='rootedTre_orthologymanyToMany.txt'

for file in d.keys(): 
    commande='python3 Nametunicates_hier.py -f OrthologyRelationships/'+d[file]+' -out Outputs/'+file+' -check n'
    os.system(commande)
    commande='python3 Nametunicates_all.py -f OrthologyRelationships/'+d[file]+' -out Outputs/'+file+' -check n'
    os.system(commande)