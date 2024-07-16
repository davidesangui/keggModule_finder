import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("annotationFolder",help="Folder containing annotations, one per genome.")
parser.add_argument("outputPrefix",help="Prefix of the output files.")
parser.add_argument('-format',choices=['list','emapper'],required=True,help='Format of the annotation files. LIST: a text file where each line lists a Kegg orthologue if the genome can encode it. EMAPPER: output of eggnog mapper, in particular the < .emapper.annotations > file. File names must end with this suffix.')
group=parser.add_mutually_exclusive_group(required=True)
group.add_argument('-definition_file',help='Text file where each line lists a module and its definition, tab-separated.')
group.add_argument('--get_definitions',help='Set this flag if the definition file has to be generated from the Kegg database. Requires Bio.KEGG package and an internet connection. The definition file is saved in < *outputPrefix*_definitions.tsv >',action='store_true')


args=parser.parse_args()

def module_length(expr):
    a=[]
    toapp=''
    opn=0
    for i in range(len(expr)): # firstly it divides in a list the different step of the module
        toapp+=expr[i]
        if expr[i]=='(':
            opn+=1
        elif expr[i]==')':
            opn-=1
        elif expr[i]==' ' and opn==0:
            a.append(toapp.strip())
            toapp=''
        if i==len(expr)-1:
            a.append(toapp)
    return(len(a))

def module_solver(expr,annotation_set): # expr is the definition of the module
    a=[]
    toapp=''
    opn=0
    for i in range(len(expr)): # firstly it divides in a list the different step of the module
        toapp+=expr[i]
        if expr[i]=='(':
            opn+=1
        elif expr[i]==')':
            opn-=1
        elif expr[i]==' ' and opn==0:
            a.append(toapp.strip())
            toapp=''
        if i==len(expr)-1:
            a.append(toapp)
    
    converted=[] # then each step is converted to a boolean string using the dictionary
    for step in a:
        conv=''
        inminus=False
        for i in range(len(step)):
            if not inminus:
                if step=='--':
                    continue
                if step[i]==' ':
                    conv+='and '
                elif step[i]=='+':
                    conv+='and '
                elif step[i]=='-' and step[i+1]=='(':
                    inminus=True
                elif step[i]==',':
                    conv+='or '
                elif step[i]=='(' or step[i]==')':
                    conv+=step[i]
                elif step[i]=='K' and step[i-1]!='-':
                    if step[i:i+6] in annotation_set:
                        conv+='True '
                    else:
                        conv+='False '                    
            else:
                if step[i]==')':
                    inminus=False
        converted.append(conv)
    
    pres=0
    for i in range(len(converted)): # then each boolean step is evaluated and the number of present step is returned
        if converted[i]=='': 
            continue
        if eval(converted[i]):
            pres+=1
    return(pres)

def module_definer(module_file):
    module_definition={} 
    f=open(module_file) 
    for line in f:
        a=line.strip().split('\t')
        module_definition[a[0]]=a[1]
    f.close()
    return(module_definition)

def annot_emapper(annot_folder):
    genome_koset={}
    files=[annot_folder+'/'+x for x in os.listdir(annot_folder) if x.endswith('.emapper.annotations')]  
    if not files: print('No eggnog output files were found in %s.\nFile names must end with < .emapper.annotations >.\nAborting...'%annot_folder); exit(1)  
    for fi in files:
        f=open(fi)
        genome='.'.join(fi.split('/')[-1].split('.')[:-2])
        genome_koset[genome]=set()
        for line in f:
            if line[0]!='#':
                fields=line.strip().split('\t')
                if fields[11]!='-': 
                    for ko in fields[11].split(','):
                        genome_koset[genome].add(ko[3:])
        f.close()
    return {x:genome_koset[x] for x in sorted(genome_koset.keys())} 

def annot_list(annot_folder):
    genome_koset={}
    files=[annot_folder+'/'+x for x in os.listdir(annot_folder)]  
    for fi in files:
        f=open(fi)
        genome='.'.join(fi.split('/')[-1].split('.')[:-1]) if len(fi.split('/')[-1].split('.'))>1 else fi.split('/')[-1]
        genome_koset[genome]={x.strip() for x in f}
        f.close()
    return {x:genome_koset[x] for x in sorted(genome_koset.keys())} 

####### main #######

if not args.get_definitions:
    definition_file=args.definition_file
else:
    from Bio import SeqIO
    from Bio.KEGG import REST
    from Bio.KEGG.KGML import KGML_parser
    print('Getting module definitions, be patient... (next time you can recycle <%s_definitions.tsv> )'%args.outputPrefix)
    listofmodules=[x.split('\t')[0] for x in REST.kegg_list('module').readlines()]
    new=open(args.outputPrefix+'_definitions.tsv','w')
    howmany=len(listofmodules)
    for s in range(howmany):
        a=REST.kegg_get(listofmodules[s]).read().split('\n')
        for fie in a:
            b=fie.split()
            if b[0]=='DEFINITION':
                defin=' '.join(b[1:])
                print(listofmodules[s]+'\t'+defin,file=new)
                break
        print('Done:',str(s+1)+'/'+str(len(listofmodules)),end='\r') if s+1<howmany else print('Done:',str(s+1)+'/'+str(len(listofmodules)))
    new.close()
    definition_file=args.outputPrefix+'_definitions.tsv'

if args.format=='list': annot_dict=annot_list(args.annotationFolder)
elif args.format=='emapper': annot_dict=annot_emapper(args.annotationFolder) 

mod_def=module_definer(definition_file)

new=open(args.outputPrefix+'_count.tsv','w')
new2=open(args.outputPrefix+'_minusOne.tsv','w')

topr='genome'
topr2='genome'
for mod in mod_def:
    topr+='\t'+mod+',length='+str(module_length(mod_def[mod]))
    topr2+='\t'+mod
print(topr,file=new)
print(topr2,file=new2)

for gen in annot_dict:
    annot_set=annot_dict[gen]
    topr=gen
    topr2=gen
    for module in mod_def:
        definition=mod_def[module]
        if 'M' in definition: # skip nested modules
            topr+='\t'+'0'
            topr2+='\t'+'0'
            continue
        ko_pres=module_solver(definition,annot_dict[gen])
        topr+='\t'+str(ko_pres)
        topr2=topr2+'\t'+'1' if ko_pres>=module_length(mod_def[mod])-1 else topr2+'\t'+'0'
    print(topr,file=new)
    print(topr2,file=new2)
new.close()
new2.close()
