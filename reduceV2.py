import argparse  
import subprocess
import os 

#Reading arguments
parser=argparse.ArgumentParser()
parser.add_argument("-f","--file",help="OrthologyRelationships file ",type=str,required=True)
args=parser.parse_args()
file=args.file

#Define the hierarchy order
Ordervert=['Hsap','Mmus','Psin','Ggal','Lcha','Cmil']

#fl is the number of same forst letter needed to be considered as the same name 
fl=3

for f in os.listdir('Named/'): 
    os.remove('Named/'+f)

#Fuction to find and return occurence of a string 
def grep(str,file): 
    commande='grep '+str+' '+file
    #print(commande)
    res=subprocess.check_output(commande.split(' '))
    return res.decode().split("\n")[:-1]   #Last element is a blanc space so lets remove it 


#Apply the hierarchy to the file
def hierV2(file): 
    with open(file) as f, open('Outputs/hierV2_'+file.split('/')[len(file.split('/'))-1],'w') as hier: 
        #initialize
        lines=f.readlines()
        current=''
        relist=[]
        for line in lines: 
            #For each line we look a the tunicate genes, if it's a new one, we look at the relationships of the previous one
            if line.split(' ')[0]!=current : 
                current=line.split(' ')[0]
                if len(relist)>=1 and relist[0]!=' ': 
                    k=0 
                    stop=0
                    #We start by searching relationship with the first vertebrate and repeat the search with the second and so one , till reach the and of the vertebrate list
                    while k<len(Ordervert) and stop==0: 
                        #print(relist)
                        for rel in relist: 
                            if rel.split(' ')[1].split('|')[0]==Ordervert[k]: 
                                stop=1
                                #print(rel)
                                #if a relationship is find with the current vertebrate we write it and set stop to 1 so we don't search relationships with other vertebrate
                                hier.write(rel)
                        if k==0 and stop==1: 
                             for rel in relist : 
                                if rel.split(' ')[1].split('|')[0]==Ordervert[1]: 
                                    hier.write(rel)
                        k+=1
                relist=[]
            #list relationship is add at each line and reset at each different tunicate genes found
            relist.append(line)
        #Don't forget to look at the final tunicate genes 
        if len(relist)>=1 and relist[0]!=' ': 
                    k=0 
                    stop=0
                    while k<len(Ordervert) and stop==0: 
                        #print(relist)
                        for rel in relist: 
                            if rel.split(' ')[1].split('|')[0]==Ordervert[k]: 
                                #print(rel)
                                stop=1
                                hier.write(rel)
                        k+=1

#Function for test if a string got a similar string (return True if the string is not a name to avoid make sure to add only names in the tab in the next function)
def intab(str,tab): 
    #Check if the string is a good name, if not return True to avoid adding this into the list 
    if str=='' or str.strip()=='notFound': 
        return True
    #Check if there is a similar string (same first letters) in the tab 
    for t in tab : 
        str=str.upper()
        t=t.upper()
        if str[:fl]==t[:fl]: 
            return True
    #If not return False :'the string is a name and as no similar string in the list ' 
    return False 

#Isolate tunicate genes where we find more than one name and genes where no name is found, and name the other ones 
def get_difname(file): 
    with open(file) as f,open('Outputs/Diffnames_'+file.split('/')[len(file.split('/'))-1],'w') as dif,open('Named/1name_'+file.split('/')[len(file.split('/'))-1],'w') as one,open('Unnamed.txt','w') as un:
        #Initialize
        lines=f.readlines()
        current=''
        current_name=[]
        current_fullname=[]
        difnames={}
        for line in lines :
            ascg=line.split(' ')[0] 
            #For each line we look a the tunicate genes, if it's a new one, we look at the relationships of the previous one
            if ascg!=current: 
                #Update the difnames dictionnary 
                if len(current_name) in difnames.keys(): 
                    difnames[len(current_name)]+=1
                else : 
                    difnames[len(current_name)]=1 
                #If we found more than one name we write in the corresponding file all the relationships of the tunicate genes
                if len(current_name)>1: 
                    gr=grep(current,file)
                    for g in gr: 
                        dif.write(g.replace('\n','')+"\n")
                #If we found just one name we name the genes
                if len(current_name)==1: 
                    one.write(current+' NameFound :'+nametab(list(set(current_fullname)))+'\n')
                #If gene can't be named we also put it in the right file 
                if len(current_name)==0 and current!='': 
                    un.write(current+'\n')
                current=ascg
                current_name=[] 
                current_fullname=[]
            #At each line the name list is update and it's reset when a new tunicate gene is found 
            vertg=line.split(' ')[1]  
            current_fullname.append(vertg.split('Name::')[1].split(';')[0])
            if not(intab(vertg.split('Name::')[1].split(';')[0],current_name)) : 
                current_name.append(vertg.split('Name::')[1].split(';')[0])
        #And don't forget the last tunicate gene 
        if len(current_name) in difnames.keys(): 
            difnames[len(current_name)]+=1
        else : 
            difnames[len(current_name)]=1 
        if len(current_name)>1: 
            gr=grep(current,file)
            for g in gr: 
                dif.write(g.replace('\n','')+"\n")
        if len(current_name)==1: 
                one.write(current+' NameFound :'+nametab(list(set(current_fullname)))+'\n')
        if len(current_name)==0: 
                un.write(current+'\n')
    f.close()
    return difnames

#Fuction to to search a similar string in a list 
def intabis(str,tab): 
    #first remove string that are not a name 
    try : 
        tab.remove('')
    except ValueError : 
        a=1
    try :
        tab.remove('None\n')
    except ValueError :
        a=1
    #Then search for a similar string (same 2 first letters) in the list
    for t in tab : 
        str=str.upper().strip()
        t=t.upper().strip() 
        if str[:fl]==t[:fl]: 
            return t
     #if a similar strinf is found, return it , else return None 
    return None

#Function to build a syndic with a name and the similar names/synonyms of the other ortologous
def findsyn(namedic):
     simidic={}
     nb=len(namedic.keys())
     for name in namedic: 
          simidic[name]=[]
          #For each name check all the name/synonyms similar in the orther ortologous 
          for n2 in namedic: 
               if intabis(name,n2.split('-'))!=None :
                    simidic[name].append(intabis(name,n2.split('-')))
                    #print(simidic)
               elif intabis(name,namedic[n2])!=None :
                    simidic[name].append(intabis(name,namedic[n2]))
          simidic[name]=list(set(simidic[name]))
     nb=nb-max([len(simidic[n]) for n in simidic]+[0])
     #nb is the minimum number of gene without a similar name/synonym so 
     # if nb=0 it means that there is a name that can linked all the genes together  
     return simidic,nb

#Same fuction but synonyms are checked at the first step 
def findsynbis(namedic):
     simidic={}
     nb=len(namedic.keys())
     for name in namedic: 
          for syn in namedic[name]: 
            simidic[syn]=[]
            for n2 in namedic: 
                if intabis(syn,n2.split('-'))!=None :
                        simidic[syn].append(intabis(syn,n2.split('-')))
                        #print(simidic)
                elif intabis(syn,namedic[n2])!=None :
                        simidic[syn].append(intabis(syn,namedic[n2]))
            simidic[syn]=list(set(simidic[syn]))
     nb=nb-max([len(simidic[n]) for n in simidic]+[0])
     return simidic,nb

#Transform a string into a tab where follwing digit are together to avoid reduction of 21 and 23 into 2 1/3
def strtotab(s): 
    t=[]
    digit=''
    #if string is only a car do nothing 
    if len(s)==1 : 
        t=[s]
    else : 
        for i in range(0,len(s)-1):  # else add digit together once the following char is not a digit 
            if s[i].isdigit(): 
                digit+=s[i]
                if not(s[i+1].isdigit()) :
                    t.append(digit)
                    digit=''
            else :                   # add char when it's not a digit
                t.append(s[i])
        # don't forget to add the last one
        if s[len(s)-1].isdigit(): 
            digit+=s[len(s)-1]
            t.append(digit)
        else : 
            t.append(s[len(s)-1])
    return t

#Transform a list of different name but similar into one name 
def nametab(tab): 
    #remove tab elements that are not a name 
    try : 
        tab.remove('')
    except ValueError : 
        a=1
    try :
        tab.remove('notFound\n')
    except ValueError : 
        a=1
    i=0
    stop=0 
    #turn string into a tab where following digit are together
    tab=[strtotab(i) for i in tab]
    # Find the commun radix by looking at the first 2 letters then the third and so one till find differences 
    while stop!=1 and i<min([len(i)for i in tab]):
        i+=1
        if i>=min([len(i)for i in tab]): 
            stop=1
        else : 
            c=tab[0][i]
            for o in tab[1:] :
                if o[i]!=c :
                    stop=1
    tabbis=[''.join(t[i:]) for t in tab]
    a=0
    while a==0 : 
        try : 
            tabbis.remove('')
        except ValueError : 
            a=1
    #Return the commun radix +/suffix for all the suffix found 
    return ''.join(tab[0][:i])+' '+'/'.join(tabbis)*(len(tab)>1)

# Look at all genes to try to find commun similar name/syn to give a name that refere all the relationships
def get_syn(file): 
    with open(file) as f,open('Outputs/Nosyn_'+file.split('/')[len(file.split('/'))-1],'w') as dif,open('Named/Withsyn_'+file.split('/')[len(file.split('/'))-1],'w') as ws:
        #Initialize
        lines=f.readlines()
        current=''
        current_name=[]
        nmdic={}
        for line in lines :
            ascg=line.split(' ')[0] 
            vertg=line.split(' ')[1]  
            #For each line we look a the tunicate genes, if it's a new one, we look at the relationships of the previous one
            if ascg!=current: 
                for lin in current_name : 
                    #Build the dictionnary where a gene name is associted with synonmys from the sam genes 
                    try : 
                        name=lin.split('Name::')[1].split(';')[0]
                        syn=lin.split('_Syn.')[1].split('::')
                        nmdic[name]=syn
                    except IndexError : 
                        a=1
                simi,nb=findsyn(nmdic)
                if nb>=1: 
                    simi,nb=findsynbis(nmdic)
                #If nb>=1 it means that there is no commun name/synonym  so we write all the relation of the tunicate gene in the right file 
                if nb>=1: 
                    for lin in current_name : 
                        dif.write(lin.replace("\n",'')+"\n")
                elif nb<1 and current!='':
                #If nb<1 it mean that there is a commun name/synonym so we name the tuncate gene using the nametab function 
                    a=0
                    simax=[]
                    for i in simi : 
                        if len(simi[i])>=a: 
                            simax=simi[i]
                            a=len(simi[i])
                    ws.write(current+' NameFound :'+nametab(simax)+'\n')
                nmdic={}
                current=ascg
                current_name=[]
            current_name.append(current+' '+vertg)
        #Don't forget the last one here again  
        for lin in current_name : 
                try : 
                    name=lin.split('Name::')[1].split(';')[0]
                    syn=lin.split('_Syn.')[1].split('::')
                    nmdic[name]=syn
                except IndexError : 
                    a=1
        simi,nb=findsyn(nmdic)
        if nb>=1: 
            simi,nb=findsynbis(nmdic)
        if nb>=1: 
            for lin in current_name : 
                dif.write(lin.replace("\n",'')+"\n")
        elif nb<1 and current!='':
            a=0
            simax=[]
            for i in simi : 
                if len(simi[i])>a: 
                    simax=simi[i]
                    a=len(simi[i])
            ws.write(current+' NameFound :'+nametab(simax)+'\n')
    f.close()

#Isolate the tunciates genes where we find more than 2 differents name (2 first letters differents) and name those where we find only 2 names
def get_more2(file): 
    with open(file) as f,open('Outputs/More2_'+file.split('/')[len(file.split('/'))-1],'w') as dif,open('Named/only2_'+file.split('/')[len(file.split('/'))-1],'w') as two:
        #Initialize
        lines=f.readlines()
        current=''
        current_name=[]
        difnames={}
        for line in lines :
            ascg=line.split(' ')[0] 
            #For each line we look a the tunicate genes, if it's a new one, we look at the relationships of the previous one
            if ascg!=current: 
                #Update the difnames dictionnary 
                if len(current_name) in difnames.keys(): 
                    difnames[len(current_name)]+=1
                else : 
                    difnames[len(current_name)]=1 
                #If we found more than 2 names we write in the corresponding file all the relationships of the tunicate genes
                if len(current_name)>2: 
                    gr=grep(current,file)
                    for g in gr: 
                        dif.write(g.replace('\n','')+"\n")
                #if we found only 2 genes we name the gene by joining the 2 names with a '-'
                if len(current_name)==2 :
                    two.write(current+' NameFound :'+current_name[0]+'-'+current_name[1]+'\n')
                current=ascg
                current_name=[] 
            #At each line the name list is update and it's reset when a new tunicate gene is found    
            vertg=line.split(' ')[1]  
            if not(intab(vertg.split('Name::')[1].split(';')[0],current_name)) : 
                current_name.append(vertg.split('Name::')[1].split(';')[0])
        #And don't forget the last tunicate gene       
        if len(current_name) in difnames.keys(): 
            difnames[len(current_name)]+=1
        else : 
            difnames[len(current_name)]=1 
        if len(current_name)>2: 
            gr=grep(current,file)
            for g in gr: 
                dif.write(g.replace('\n','')+"\n")
        if len(current_name)==2 :
                    two.write(current+' NameFound :'+current_name[0]+'-'+current_name[1]+'\n')
    f.close()
    return difnames

#Function for check if two list countains the same elements 
def same_elem(tab1,tab2): 
    a=[i for i in tab2]
    for i in tab1: 
        #for each element in the first list we check if it is in the second 
        if not(i in tab2): 
            return False 
        else: 
            #we remove it from the copy list of the second list to make sure that all elements are in the same number between the 2 lists
            try : 
                a.remove(i)
            except ValueError : 
                return False
    #If each elements of the first list are in the second at the same number and the second list contains no other elements return true 
    if len(a)==0 : 
        return True
    else : 
        return False

#Check if a list of list contain a list with same elements than a list 
def tintt(t,tt): 
    #Make sure to not add an empty list in the tab
    if t==[]: 
        return True 
    for i in tt: 
        if same_elem(t,i): 
            return True 
    return False 

#Compute the numbe of differente case(=Same ortologous)
def get_nb_case(file): 
    with open(file) as f,open('Outputs/diffCase_'+file.split('/')[len(file.split('/'))-1],'w') as dif,open('Tocomplex.txt','w') as tc:
        #Initialize
        lines=f.readlines()
        current=''
        current_name=[]
        cases=[]
        for line in lines :
            ascg=line.split(' ')[0] 
            vertg=line.split(' ')[1]
            #For each line we look a the tunicate genes, if it's a new one, we look at the relationships of the previous one
            if ascg!=current : 
                tc.write(ascg)
                #Check if the naming case have already been seen or not 
                if not(tintt(current_name,cases)) : 
                    cases.append(current_name)
                    #If it's a new case write the first example of the cases in the right file 
                    for g in current_name: 
                        dif.write(current+' '+g.replace('\n','')+'\n')
                current=ascg
                current_name=[]
            #At each line the name list is update and it's reset when a new tunicate gene is found 
            current_name.append(vertg)
        #And don't forget the last tunicate gene 
        if not(tintt(current_name,cases)): 
            cases.append(current_name)
            for g in current_name: 
                dif.write(current+' '+g.replace('\n','')+'\n')
        return len(cases)

#Apply the function wrote before 
print('Proceding...')
hierV2(file)
print('hier done')
get_difname('Outputs/hierV2_'+file.split('/')[len(file.split('/'))-1])
print('difname done ')
get_syn('Outputs/Diffnames_hierV2_'+file.split('/')[len(file.split('/'))-1])
print('syn done ')
get_more2('Outputs/Nosyn_Diffnames_hierV2_'+file.split('/')[len(file.split('/'))-1])
print('more2 done ')
print('Number of cases left : '+str(get_nb_case('Outputs/More2_Nosyn_Diffnames_hierV2_'+file.split('/')[len(file.split('/'))-1])))
print('End')