import subprocess
import argparse 
import matplotlib.pyplot as plt 

#Fuction to find and return occurence of a string 
def grep(str,file): 
    commande='grep '+str+' '+file
    #print(commande)
    res=subprocess.check_output(commande.split(' '))
    return res.decode().split("\n")[:-1]   # dernière element à enlever car toujours un espaces blanc 



#Reading argument 
parser=argparse.ArgumentParser()
parser.add_argument("-f","--file",help="OrthologyRelationships file ",type=str,required=True)
parser.add_argument("-show","--Stats",help="Stats to compute",type=str,required=False,default='all',choices=
                    ['all','Nbtun','Nbtunbysp','Nbrel','Nbuniquerel','Nbnf','Nbdiffname','Nbcase'])
parser.add_argument("-gr","--graph",help="Choose if you want to plot graphics or not",required=False,default='Y',choices=['Y','n'])

args=parser.parse_args()
file=args.file
stats=args.Stats.split(' ')
graph=args.graph



def draw_pie(results,names,title): 
    plt.figure(figsize=(5,5))
    plt.pie(results, labels = names,
           autopct = lambda x: str(round(x, 1)) + '%',
           pctdistance = 0.7, labeldistance = 1.1,
           shadow = True)
    plt.title(title)
    #plt.show()

#Compute the number of tunicate gene models in the file 
def get_nb_asci(file):
    with open(file) as f :
        lines=f.readlines()
        current=''
        n=0
        for line in lines :
            ascg=line.split(' ')[0]
            if current!=ascg:
                current=ascg
                n+=1
    
    f.close()
    return n

#Compute the number of tunicate gene models in the file by species
def get_ascistats(file): 
    asc=['Phmamm','Phfumi','Cisavi','Cirobu','Moocci','Moocul','Mooccu','Boschl','Boleac','Haaura','Harore','Coinfl','Stclav']
    nbasc={}
    for a in asc : 
        nbasc[a]=0
    with open(file) as f :
        lines=f.readlines()
        current=''
        for line in lines :
            ascg=line.split(' ')[0]
            if current!=ascg:
                current=ascg
                nbasc[ascg.split('|')[0]]+=1
    return nbasc


#Compute the nomber of relation by vertebrates
def get_nbrel(file):
    Ordervert=['Hsap','Mmus','Psin','Ggal','Lcha','Cmil']
    with open(file) as f :
        lines=f.readlines()
        current=''
        dejavu={}
        tot={}
        for v in Ordervert :  
            dejavu[v]=0
            tot[v]=0
        for line in lines :
            ascg=line.split(' ')[0]
            vertg=line.split(' ')[1]
            if current!=ascg:
                current=ascg
                #If it's a new tunicate gene models add all his vertebrate relationships to the total
                for v in Ordervert : 
                    tot[v]+=dejavu[v]
                    dejavu[v]=0
            for v in Ordervert :
                if vertg.split('|')[0]==v: 
                    dejavu[v]+=1 
        for v in Ordervert :    # Don't forget the last one 
                tot[v]+=dejavu[v]
                dejavu[v]=0
    f.close()
    return tot

#Compute the number of unique relationships by vertebrate (no more than 1 per tunicate gens )
def get_nbunique_rel(file):
    Ordervert=['Hsap','Mmus','Psin','Ggal','Lcha','Cmil']
    with open(file) as f :
        lines=f.readlines()
        current=''
        dejavu={}
        tot={}
        for v in Ordervert :  
            dejavu[v]=0
            tot[v]=0
        for line in lines : 
            ascg=line.split(' ')[0]
            vertg=line.split(' ')[1]
            if current!=ascg:
                current=ascg
                #If it's a new tunicate gene models add all his vertebrate relationships to the total
                for v in Ordervert : 
                    tot[v]+=dejavu[v]
                    dejavu[v]=0
            for v in Ordervert  :  
                if vertg.split('|')[0]==v: 
                    dejavu[v]=1 
        for v in Ordervert : # Don't forget the last one 
            tot[v]+=dejavu[v]
    f.close()
    return tot

#Compute the number of notfound and no name in the file 
def get_nf_and_noname(file):
    Ordervert=['Hsap','Mmus','Psin','Ggal','Lcha','Cmil']
    with open(file) as f :
        lines=f.readlines()
        tot={}
        nf={}
        noname={}
        for v in Ordervert :  
            tot[v]=[]
            noname[v]=[]
            nf[v]=[]
        #Look at each line and add +1 in nf if this is a notfound, +1 in noname if name is empty 
        for line in lines : 
            vertg=line.split(' ')[1]
            for v in Ordervert : 
                if vertg.split('|')[0]==v : 
                    if vertg not in tot[v] : 
                        tot[v].append(vertg)
                        if vertg.split('Name::')[1].split(';')[0]=='': 
                            noname[v].append(vertg)
                        if vertg.split('Name::')[1].split(';')[0].strip()=='notFound': 
                            nf[v].append(vertg)     
    f.close()
    print("Number of vertebrate genes : ")
    for v in Ordervert : 
        print(f"{v} :  Total = {len(tot[v])}  Notfound = {len(nf[v])}  Unnamed = {len(noname[v])}")
    return tot,nf,noname


def draw_bar(results,names,title): 
        #x=[(len(nf[v])+len(noname[v]))/len(tot[v]) for v in tot]
    plt.figure(figsize=(5,5))
    plt.bar([i for i in range(len(results))],results,width=0.5,color='indigo')
    plt.xticks([i for i in range(len(results))],names)
    plt.title(title)

#Function to test if a string is in a tab or not 
def intab(str,tab): 
    if str=='' or str.strip()=='notFound': 
        return True
    for t in tab : 
        str=str.upper()
        t=t.upper()
        #if str==t or str[:-1]==t or str==t[:-1] or str[:-1]==t[:-1]: 
        if str[:2]==t[:2]: 
            return True
    #print(str,tab)    
    return False 

#Compute the distribution of the number of differentes name found in vertebrate ortologous of a same tunicate gene 
def get_difname(file): 
    with open(file) as f:
        lines=f.readlines()
        current=''
        current_name=[]
        difnames={}
        for line in lines :
            ascg=line.split(' ')[0] 
            #if the number of differents names is tottally new create the key corresponding if not just add +1
            if ascg!=current: 
                if len(current_name) in difnames.keys(): 
                    difnames[len(current_name)]+=1
                else : 
                    difnames[len(current_name)]=1 
                current=ascg
                current_name=[] 
            vertg=line.split(' ')[1] 
            #Add the new name at each lines (if it's a new one) 
            if not(intab(vertg.split('Name::')[1].split(';')[0],current_name)) : 
                current_name.append(vertg.split('Name::')[1].split(';')[0])
        #Don't forget the last one 
        if len(current_name) in difnames.keys(): 
            difnames[len(current_name)]+=1
        else : 
            difnames[len(current_name)]=1 
    f.close()
    return difnames

#Fuction to check if two list got the same elements 
def same_elem(tab1,tab2): 
    a=[i for i in tab2]
    for i in tab1: 
        if not(i in tab2): 
            return False 
        else: 
            try : 
                a.remove(i)
            except ValueError : 
                return False
    if len(a)==0 : 
        return True
    else : 
        return False


#fuction to check if a list is in a list of list 
def tintt(t,tt): 
    if t==[]: 
        return True
    for i in tt: 
        if same_elem(t,i): 
            return True 
    return False 

#Compute the number of case (=naming problem or orthologous group) int the file 
def get_nb_case(file): 
    with open(file) as f:
        lines=f.readlines()
        current=''
        current_name=[]
        cases=[]
        for line in lines :
            ascg=line.split(' ')[0] 
            vertg=line.split(' ')[1]
            if ascg!=current : 
                #print(current_name)
                #For eahc new tunicate genes, check the naming case , if this is a new one +1 and add to cases 
                if not(tintt(current_name,cases)) : 
                    cases.append(current_name)
                current=ascg
                current_name=[]
            current_name.append(vertg)
        if not(tintt(current_name,cases)): 
            cases.append(current_name)
        return len(cases) 


#Apply the function 

if 'Nbtun' in stats or 'all' in stats : 
    tot=get_nb_asci(file)
    print("Number of tunicates genes : "+str(tot)+"\n")

if 'Nbtunbysp' in stats or 'all' in stats : 
    tot=get_nb_asci(file)
    dic=get_ascistats(file)
    print("Number of tunicate genes by species : "+str(dic)+"\n")
    if graph=='Y': 
        draw_pie([dic[t]/tot for t in dic],dic.keys(),'Distribution of gene number by tunicate species')

if 'Nbrel' in stats or 'all' in stats : 
    dic=get_nbrel(file)
    print("Number of Relationships by vertebrate : "+str(dic)+"\n")
    if graph=='Y': 
        tot=sum([dic[v] for v in dic])
        draw_pie([dic[v]/tot for v in dic],dic.keys(),'Distribution of relationships by vertebrates')


if 'Nbuniquerel' in stats or 'all' in stats : 
    print("Nombre of tunicate gene with at least one relation with this vertebrate : "+str(get_nbunique_rel(file))+'\n')

if 'Nbnf' in stats or 'all' in stats : 
    tot,nf,noname=get_nf_and_noname(file)
    print("\n")
    if graph=='Y' : 
        draw_bar([(len(nf[v])+len(noname[v]))*100/len(tot[v]) for v in tot ],[v for v in tot],'Percentage of Unfound or Unnamed genes by vertebrate')

if 'Nbdiffname' in stats or 'all' in stats: 
    nbdifname=get_difname(file)
    if graph=='Y': 
        if len(nbdifname.keys())>3: 
            ndif=[nbdifname[i] for i in range(0,4)]
            for i in nbdifname:
                if i>3: 
                    ndif[3]+=nbdifname[i]
            l=[str(i) for i in range(0,3)]
            l.append('3 or more')
        else : 
            l=nbdifname.keys()
            ndif=[nbdifname[i] for i in nbdifname]
        draw_pie(ndif,l,'Distribution of Number of differents name')

if 'Nbcase' in stats or 'all' in stats : 
    print("Number of different naming case : "+str(get_nb_case(file)))

plt.show()

