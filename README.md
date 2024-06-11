
Context : 

Being part of a tunicates studiyng team, this work aim to give tunicates genes a symbol using orthologies relationships with vertebrates. 
Previous studies had etablished those relationships with six vertebrates, Human, Mouse,Chicken,Chinese Turtle,coelcanth and Australian ghostshrak. 
All tunicates genes got a Gene model indentifier etablished using the rules describe in Guidelines for the Nomenclature of Gentetic Elements in Tuncates Genomes (Lemaire et al.2015). 
This identifier gave information about the assembly, the location and the instance of the gene. 
Give those genes a primary symbol also respecting the guidelines will be usefull to have information about the gene function and ortohlogie relationships to make futher tunicates studies much easier.


Methodology : 

Give a symbol to gene using orthology relationships is a really specific problem and as there a lot of genes here the goal is to make this ‘naming’ automatic. 
In a very simple way this work use ortohology relationships file and out file where every tunicates genes with at least one relations with a vertebrate is associated to a symbol. 
But some gene are harder to ‘name’ because of the variety of symbol founds in the vertebrates ortologous. Method used here start by etablished a hierarchy between vertebrates to reduced the variety of name.
Human’s orthologous are used at first, if a gene doesn’t have one then mouse ones are used and so on respecting the following line: 

Human→Mouse→Chinese turtle → Chicken → Coelcanth → Australian ghostshark. 

But there still gene hard to ‘name. So once gene with really different name (not same two first letters), the third step, still following the guidelines of the article, is to try to find the Ciona intestatinalis ortologous to use it symbol. 
For genes who havn’t got one with a symbol established, a  commun name or radical is searched between first symbols and synonyms of each gene in the ortology relationship. 
Here a few genes left and some juste have 2 different symbol so a fusion of those two names is used but when there is more than 2, gene symbols are established manually doing some research to find the best name that will refer to all ortologous, for example used the name of the family whenever this symbol is not so large and there is almost all the familiy.


Tools and implementation : 

All the work done in this project had been done using python (version 3.10.12)  in visual studio code. The main reason being that I very familiar to this environnement and parsing text files with python is simple. 

Input data must be text file with relationships wrote in the following way (take care about the white space) : 

geneSpecies|geneModelID ortologousSpecies|ortologousgeneModelID

First the file Nametunicates_all.py by enter this commande : 

python3 Nametunicates_all.py -f [filename] -outf [outfilename]

If needed it is possible to only look at the first x tunicates by add -lim x. 
Terminal will show an estimated time based on previous execution time, and ask if it’s ok to continu or not. This check can be disabled by add -check n. 

The outfile contain all the ortology relationships bewteen ascidies and vertebrates in the input file and his name [outfilename]_all.txt

Then to name all the tunicate gene use the script reduce.py in the same way : 

python3 reduce.py -f [filename] 

Make sure to have a folder named ‘Named’ once the script has finished the folder will countains all the tunicate genes with the name found, one file for each step of naming. 

To get statistics in a file the script Stats.py can also be executed : 

python3 Stats.py -f [filename]

Stats can be configured (see the arguements help and graph of Stats) 

To add a species of interest there is three list at the beginning of the Nametunicates_all.py file (line 17-18-19)
just add the species ID as it's found in your file in the right list. 
