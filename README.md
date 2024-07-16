
Context : 

Being part of a tunicates studiyng team, this work aim to give tunicates genes a symbol using orthologies relationships with vertebrates. Previous studies had etablished those relationships with six vertebrates, Human, Mouse,Chicken,Chinese Turtle,coelcanth and Australian ghostshrak. All tunicates genes got a Gene model indentifier etablished using the rules describe in Guidelines for the Nomenclature of Gentetic Elements in Tuncates Genomes (Lemaire et al.2015). This identifier gave information about the assembly, the location and the instance of the gene. Give those genes a primary symbol also respecting the guidelines will be usefull to have information about the gene function and ortohlogie relationships to make futher tunicates studies much easier.


Methodology : 

Give a symbol to gene using orthology relationships is a really specific problem and as there a lot of genes here the goal is to make this ‘naming’ automatic. In a very simple way this work use ortohology relationships file and out file where every tunicates genes with at least one relations with a vertebrate is associated to a symbol. But some gene are harder to ‘name’ because of the variety of symbol founds in the vertebrates ortologous. Method used here start by etablished a hierarchy between vertebrates to reduced the variety of name. Human and  mouse orthologous are used at first, if a gene doesn’t have one then the chinese turtle ones are used and so on respecting the following line: 

Human&Mouse→Chinese turtle → Chicken → Coelcanth → Australian ghostshark. 

But there still gene hard to ‘name’. For those genes a commun name or radical is searched between first symbols and synonyms of each gene in the ortology relationship. Here a few genes left and some juste have 2 different symbol so a fusion of those two names is used but when there is more than 2, gene symbols are established manually doing some research to find the best name that will refer to all ortologous, for example used the name of the family whenever this symbol is not so large and there is almost all the familiy.


Tools and implementation : 

All the work done in this project had been done using python (version 3.10.12)  in visual studio code. 

Input data must be text file with relationships wrote in the following way (take care about the white space) : 

geneSpecies|geneModelID ortologousSpecies|ortologousgeneModelID

All relationships must be in the same directory than countains only those relationships ([Relationships directory]). This directory must be in a directory with the 5 script files ‘Nameall.py’, ’Nametunicates_all.py’, ’reduceV2.py’, ‘Correctname.py’ and ‘Stats.py’
and the Ensembl databases for each vertebrate named by the species abreviation +.csv exemple Hsap.csv or Mouse.csv (take care Mouse symbol need to be all in upper letter!)

Then use the following commande once placed in the big directory : 

python3 Nameall.py -d [Relationships directory]

If you just want the Finalname file and no intermediate files add -r Y.

The script will first use the Nametunicates_all.py file, if needed it is possible to only look at the first x tunicates (in the relationships file) by add -lim x in the call of Nametunicates_all.py.
Don’t worry if there is no x genes model at the , it just means that some of those genes don’t get relationships with vertebrate in the dataset so they are not showed. 

Terminal will show an estimated time based on previous execution time, and ask if it’s ok to continu or not. This check can be disabled by add -check n. 

LEGROS_all.txt countains now all the relationships tunicates-vertebrates if your dataset and 

Then the script will call reduceV2.py at the end of the execution of the file if there is some naming problem unresolved script will show a list of the name found for each cases 
write - to skip the cases , -exit to skip all the following cases or the name wanted

Once Correctname.py is executed the file Finalname countain each tunicate genes named with his name found 

To get statistics in a file the script Stats.py can also be executed : 

python3 Stats.py -f [filename]

Stats can be configured (see the arguements help of graph and Stats) 
