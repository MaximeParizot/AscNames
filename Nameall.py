import argparse 
import os 
import time 

a=time.time()
parser=argparse.ArgumentParser()
parser.add_argument("-d","--directory",help="OrthologyRelationships directory ",type=str,required=True)
parser.add_argument("-r","--remove",help=' Remove intermediate files ? Y/n',type=str,choices=['Y','n'],default='n',required=False)
parser.add_argument("-check","--checktime",help="limit the number of tunicates (used for test)",type=str,choices=['Y','n'],default='Y',required=False)
args=parser.parse_args()
check=args.checktime
dir=args.directory 
rem=args.remove

if not(os.path.exists('Outputs/')) : 
    os.mkdir('Outputs')

n=0
with open('Outputs/LEGROS.txt','w') as lg : 
    for file in os.listdir(dir) :
        with open(dir+'/'+file) as f : 
            for line in f.readlines(): 
                lg.write(line)
                n+=1


if check=='Y': 
    print(f'Estimated time in min {n*0.00028433:.3}')
    proced=''
    while proced!='Y' and proced !='n': 
        proced=input('Do you still want to proceed?[Y/n]')
    if proced=='n': 
        exit()

print('Naming vertebrate genes')
commande='python3 Nametunicates_all.py -f Outputs/LEGROS.txt -out Outputs/LEGROS -check n -lim 100 '
os.system(commande)
print('Naming of vertebrate gene completed')

if not(os.path.exists('Named/')) : 
    os.mkdir('Named')

print('reducing')
commande='python3 reduceV2.py -f Outputs/LEGROS_all.txt'
os.system(commande)
print('reducing completed')

print('correcting name')
commande='python3 Correctname.py -d Named/'
os.system(commande)
print('correcting name completed')

if rem=='Y': 
    for f in os.listdir('Named/'): 
        os.remove('Named/'+f)
    for f in os.listdir('Outputs/'): 
        os.remove('Outputs/'+f)
    os.remove('FusionNamedfile.txt')
    os.rmdir('Named')
    os.rmdir('Outputs')

print('Naming completed see Finalname.txt')
print(f'{(time.time()-a)/60:.3} minuts took')