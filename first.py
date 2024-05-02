import os

os.chdir('OrthologyRelationships')
filelist=[i for i in os.listdir()]



for file in filelist : 
    max=1
    print(file)
    with open (file) as f : 

        lines= f.readlines()
        for line in lines : 
            if max!=0: 
                print(line)
                max=max-1
    f.close()