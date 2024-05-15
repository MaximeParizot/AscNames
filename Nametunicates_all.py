import argparse 
import os 
import subprocess

parser=argparse.ArgumentParser()
parser.add_argument("-f","--file",help="OrthologieRelationships file ",type=str,required=True)
parser.add_argument("-outf","--outfile",help="Desired name for the new files ",type=str,required=True)
parser.add_argument("-lim","--limite",help="limit the number of tunicates (used for test)",type=int,default=0,required=False)
parser.add_argument("-check","--checktime",help="limit the number of tunicates (used for test)",type=str,choices=['Y','n'],default='Y',required=False)
args=parser.parse_args()
file=args.file
out=args.outfile
lim=args.limite
check=args.checktime

tun=['Phmamm','Phfumi','Cisavi','Cirobu','Moocci','Moocul','Mooccu','Boschl','Boleac','Haaura','Harore']
vert=['Cmil','Lcha','Hsap','Mmus','Ggal','Psin']
outgroupt=['Spur','Apla','Bbel','Blan']

#os.chdir('OrthologyRelationships')
#Check the file formatting
with open (file) as f : 
    try : 
        lines= f.readlines()
        for line in lines : 
            if len(line.split(' '))!=2:   # Every lines split well in 2 by ' '  
                raise Exception(ValueError,line)
            if len(line.split(' ')[0].split('|'))!=2 or len(line.split(' ')[1].split('|'))!=2 : 
                raise Exception(ValueError,line) # Every lines split well by | after first split 
    except Exception as e: 
        print(" File format Error, Line do not split well : ")
        print(e.args[1])
    f.close()

 # First build a list of tunicates genes in the file 
tunlist=[]
with open (file) as f : 
        lines= f.readlines()
        for line in lines : 
                for gene in line.split(' '): 
                    if gene.split('|')[0] in tun : 
                            tunlist.append(gene.strip())           
f.close()
tunlist=list(set(tunlist))      # To get only unique values 
#print(tunlist)
if lim!=0 :
    tunlist=tunlist[:lim]
if check=='Y': 
    #print(check)
    print(f'Estimated time in min : {len(tunlist)*10/100/60:.6f}')
    proced=''
    while proced!='Y' and proced !='n': 
        proced=input('Do you still want to proceed?[Y/n]')
    if proced=='n': 
        exit()

print('Proceding...')

#Then a dictionnary of relationships they got : 

dic={}
for g in tunlist : 
    commande='grep '+g.split('|')[1]+' '+file
    res=subprocess.check_output(commande.split(' '))
    l=res.decode().replace("\n",' ').split(' ')[:-1]   # Last element is a '' so let's remove it 
    l=list(set(l))
    l.remove(g)    # Remove also the gene himself 
    dic[g]=l
#print(dic)  


#define a fonction for more visibility : 
# take on speccies and assemblage number inarguments and return a dic with info find in Enembl
def formatout(espece,assembl): 
    #Resolving the chicken issue , dataset is different so here is a specific parser
    ##
    if espece=='Ggal': 
        commande='grep '+assembl+' Ensembl/'+espece+'.gff3'
        res=subprocess.check_output(commande.split(' ')).decode().split("\n")[:-1]
        r=[]
        for line in res : 
            if line.split("ID=")[1].split(':')[0]=='gene': 
                #print(line)
                l=line.split("ID=gene:")[1].split(';')[0]
                l=l+','+line.split("ID=gene:")[1].split(';')[0]+'.'+line.split("version=")[1].split(';')[0]
                try : 
                    l=l+','+line.split("Name=")[1].split(';')[0]+','
                except IndexError: 
                    l=l+',,'
                r.append(l)
    ## End of specific parser 
    else : 
        commande='grep '+assembl+' Ensembl/'+espece+'.csv'
        r=subprocess.check_output(commande.split(' '))
        r=r.decode().split("\n")[:-1]   # dernière element à enlever car toujours un espaces blanc 
    dic={} 
    for res in r : 
        if res.split(',')[0] in dic.keys(): 
            if res.split(',')[3]=='': 
                new=['None']
            else : 
                new=[res.split(',')[3]]

            if '.'+res.split(',')[1].split('.')[1] in sousdic[res.split(',')[2]].keys():
                 sousdic[res.split(',')[2]]['.'+res.split(',')[1].split('.')[1]].append(new[0])
            else : 
                sousdic[res.split(',')[2]]['.'+res.split(',')[1].split('.')[1]]=new
        else : 
            sousdic={}
            syn={}
            if res.split(',')[3]=='': 
                syn['.'+res.split(',')[1].split('.')[1]]=['None']
            else : 
                syn['.'+res.split(',')[1].split('.')[1]]=[res.split(',')[3]]
            sousdic[res.split(',')[2]]=syn
        dic[assembl]=sousdic
    return dic 

with open(out+'_all.txt','w') as vertf : 
    for g in tunlist : 
        res=dic[g]
        for rel in res : 
            if rel.split('|')[0] in vert:
                vertf.write(g+' '+rel)
                try : 
                    res=formatout(rel.split('|')[0],rel.split('|')[1].split('.')[0])[rel.split('|')[1].split('.')[0]]
                    for name in res.keys(): 
                          vertf.write('_Name::'+name)
                          for synid in res[name].keys(): 
                              vertf.write(';_Syn'+synid)
                              for syn in res[name][synid]: 
                                  vertf.write('::'+syn )
                    vertf.write("\n")
                except subprocess.CalledProcessError :
                    vertf.write('_Name::notFound'+"\n") 
                    print('Warning '+rel.split('|')[1].split('.')[0]+' not found in '+rel.split('|')[0])

print('End')