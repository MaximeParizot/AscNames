import subprocess
import os 
import argparse 

#Reading arguments 
parser=argparse.ArgumentParser()
parser.add_argument("-d","--dir",help="directory of Named file ",type=str,required=True)
args=parser.parse_args()
dir=args.dir

#Fuction to find and return occurence of a string 
def grep(str,file): 
    commande='grep|'+str.strip()+'|'+file
    #print(commande.split('|'))
    res=subprocess.check_output(commande.split('|'))
    return res.decode().split("\n")[:-1]   # last element to remove cause it's a white space 


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
    vert=['Hsap','Mmus','Psin','Ggal','Lcha','Cmil']
    #First find all the name with the same radical within the same species
    while vert!=[]: 
        #try the first vertebrate of the list , if not remove the first and look at the new first one 
        try : 
            namelist=list(set([i.split(',')[2] for i in grep(rad,'Ensembl/'+vert[0]+'.csv')]))
            #If name is found in human , also look at the mouse ones 
            if vert[0]=='Hsap' : 
                try :
                    namelist+=list(set([i.split(',')[2] for i in grep(rad,'Ensembl/Mmus.csv')]))
                    namelist=list(set(namelist))
                except subprocess.CalledProcessError: 
                    a=1
            vert=[]
        except subprocess.CalledProcessError :
            vert.remove(vert[0])
    #Then try to reduce the name by remove all the suffixes  found in the name 
    for end in suff:                                     
        try : 
            namelist.remove(rad+''.join(end[cpt:]))  
        except ValueError: 
            return ''.join(suff[0][:cpt])
    # if there is no suffixe left it means that the name coutain all the version of the gene so name can be reduced
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
    
#Apply the suffixe reduction to a name 
def combine(name): 
    #If there is no suffixe do nothing
    if len(name.split(' '))==1 : 
        return name[0].upper()+name[1:].lower()
    #Else create a suffixe list in the tab format following numbers are together 
    suff=name.split(' ')[1].split('/')
    stab=[]
    for s  in suff :
        stab.append(strtotab(s))
    return name.split(' ')[0][0].upper()+name.split(' ')[0][1:].lower()+reduce_suf(name.split(' ')[0],stab,0).lower()


#Function to get all the name given in a file 
def getnames(file): 
    tun=['Phmamm','Phfumi','Cisavi','Cirobu','Moocci','Moocul','Mooccu','Boschl','Boleac','Haaura','Harore','Coinfl','Stclav']
    with open(file) as f : 
        names=[]
        bigd={}
        #Look at eahc line 
        for line in f.readlines(): 
            name=line.split('NameFound :')[1].split(';')[0]
            if not(name in names): 
                #if a new name is found regroup in a dic all the tunicate genes with this name , by species 
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


#add .{a-z} and apply the reduction for all names given in eahc file of the Named directory 
def ABC(dir):
    alphabet=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    #Regroup all the named files 
    with open('Finalname.txt','w') as f , open('FusionNamedfile.txt','w') as fnf,open('Verylong.txt','w') as vl: 
        for file in os.listdir(dir): 
            with open(dir+file) as nf : 
                for line in nf.readlines(): 
                    fnf.write(line)
        fnf.close()
        #Consctruct the dictionnary of named given 
        bigd=getnames('FusionNamedfile.txt')
        for name in bigd: 
            #print(name)
            #Apply the reduction 
            fname=combine(name)
            #Check if the name is very long (more than 10 /)
            if len(fname.split('/'))>10:
                for species in bigd[name]: 
                    for i in range(0,len(bigd[name][species])) : 
                        vl.write(bigd[name][species][i])
                        vl.write('\t')
                        vl.write(fname+'\n')
            #add .{a-z} for gene model with the same name in the same species 
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
                    # just write the name found if there is only one gene model in the species with this named 
                    elif bigd[name][species]!=[]: 
                        f.write(bigd[name][species][0])
                        f.write('\t')
                        f.write(fname+'\n')

print('Proceding...')
ABC(dir)