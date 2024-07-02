import subprocess
import os 
import argparse 

parser=argparse.ArgumentParser()
parser.add_argument("-d","--dir",help="directory of Named file ",type=str,required=True)
args=parser.parse_args()
dir=args.dir

def grep(str,file): 
    commande='grep|'+str.strip()+'|'+file
    #print(commande.split('|'))
    res=subprocess.check_output(commande.split('|'))
    return res.decode().split("\n")[:-1]   # last element to remove cause it's a white space 


Ordervert=['Cmil','Lcha','Hsap','Mmus','Ggal','Psin']
tun=['Phmamm','Phfumi','Cisavi','Cirobu','Moocci','Moocul','Mooccu','Boschl','Boleac','Haaura','Harore','Coinfl','Stclav']

#Convert a string into a tab where folliwing numbers are together 
def strtotab(s): 
    t=[]
    digit=''
    if len(s)==1 : 
        t=[s]
    else : 
        for i in range(0,len(s)-1):
            if s[i].isdigit(): 
                digit+=s[i]
                if not(s[i+1].isdigit()) :
                    t.append(digit)
                    digit=''
            else : 
                t.append(s[i])
        if s[len(s)-1].isdigit(): 
            digit+=s[len(s)-1]
            t.append(digit)
        else : 
            t.append(s[len(s)-1])
    return t

#Fonction for reduce the suffixe
def reduce_suf(rad,suff,cpt): 
    vert=[i for i in Ordervert]
    #First find all the name with the same radical within the same species
    while vert!=[]: 
        try : 
            namelist=list(set([i.split(',')[2] for i in grep(rad,'Ensembl/'+vert[0]+'.csv')]))
            if vert[0]=='Hsap' : 
                try :
                    namelist+=list(set([i.split(',')[2] for i in grep(rad,'Ensembl/Mmus.csv')]))
                    namelist=list(set(namelist))
                except subprocess.CalledProcessError: 
                    a=1
            vert=[]
        except subprocess.CalledProcessError :
            vert.remove(vert[0])
    #Then 
    for end in suff:                                     
        try : 
            namelist.remove(rad+''.join(end[cpt:]))  #remove suffixes that are in the name
        except ValueError: 
            return ''.join(suff[0][:cpt])
    if len(namelist)==0: 
        return ''.join(suff[0][:cpt])
    else: 
        try : 
            suffdic={}               #else try to regroupe suffixes with the same start 
            for suf1 in suff :
                suffdic[suf1[cpt]]=[suf1] 
                for suf2 in suff : 
                    if suf1[cpt]==suf2[cpt] and suf1!=suf2: 
                        suffdic[suf1[cpt]].append(suf2)
        except IndexError : 
            return suff[0][cpt-1]
        return '/'.join([reduce_suf(rad+g,suffdic[g],cpt+1) for g in suffdic])
    

def combine(name): 
    if len(name.split(' '))==1 : 
        return name[0].upper()+name[1:].lower()
    suff=name.split(' ')[1].split('/')
    stab=[]
    for s  in suff :
        stab.append(strtotab(s))
    return name.split(' ')[0][0].upper()+name.split(' ')[0][1:].lower()+reduce_suf(name.split(' ')[0],stab,0).lower()


#Function to get all the name given in a file 
def getnames(file): 
    with open(file) as f : 
        names=[]
        bigd={}
        for line in f.readlines(): 
            name=line.split('NameFound :')[1].split(';')[0]
            if not(name in names): 
                names.append(name)
                littled={}
                for t in tun : 
                    littled[t]=[]
                for rel in grep(name,file):
                    n=rel.split('NameFound :')[1].split(';')[0]
                    if n.strip()==name.strip(): 
                        littled[rel.split('|')[0]].append(rel.split('|')[1].split('Name')[0])
                bigd[name.strip()]=littled
    return bigd



def ABC(dir):
    alphabet=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    with open('Finalname.txt','w') as f , open('FusionNamedfile.txt','w') as fnf,open('Verylong.txt','w') as vl: 
        for file in os.listdir(dir): 
            with open(dir+file) as nf : 
                for line in nf.readlines(): 
                    fnf.write(line)
        fnf.close()
        bigd=getnames('FusionNamedfile.txt')
        for name in bigd: 
            #print(name)
            fname=combine(name)
            if len(fname.split('/'))>10:
                for species in bigd[name]: 
                    for i in range(0,len(bigd[name][species])) : 
                        vl.write(bigd[name][species][i])
                        vl.write('\t')
                        vl.write(fname+'\n')
            else :
                for species in bigd[name]: 
                    if len(bigd[name][species])>=2:
                        for i in range(0,len(bigd[name][species])) : 
                            f.write(bigd[name][species][i])
                            f.write('\t')
                            if i>=len(alphabet): 
                                f.write(fname+'.'+str.lower(alphabet[int(i/len(alphabet))-1])+str(int(i-len(alphabet))%10+1)+'\n')
                            else :
                                f.write(fname+'.'+str.lower(alphabet[i])+'\n')
                    elif bigd[name][species]!=[]: 
                        f.write(bigd[name][species][0])
                        f.write('\t')
                        f.write(fname+'\n')

print('Proceding...')
ABC(dir)