## *Crenicichla tuca* genome assembly and analysis
This branch contains scripts used for all assembly steps and analysis presented in the paper presenting *C. tuca* genome.

<img src="https://github.com/mylena-s/assembly_project/blob/CtucaBranch/github.png?raw=true" width="1000">

**Reference**

All analysis were performed integrated in nextflow these workflows, except for:
* D-genies: run using an interactive local server.
* TOGA annotation: custom script runToga.sh
* FCS-NCBI: run on Usegalaxy.org
* OMark: run on https://omark.omabrowser.org/
* Small file manipulation scripts: python or bash scripts used for small taks can be found on /scripts

### Example scripts 
Examples of bash scripts used for executing the pipeline are found in /examples

All software was installed using mamba or otherwise stated in software/install.sh

# Scripts explanation
| Script | Uses |
| -------- | ------------------------------------------|
| readQC.nf | Read quality assessment before assembly |
| Assembly.nf | Genome assembly with ONT and HIFI reads |
| AssemblyCuration.nf | Scaffolding |
| AssemblyQC.nf | Assembly quality assessment, including BUSCO, kat, and other software |
| Annotation.nf | Gene and repeat annotation |	
| Centromere.nf | Centromere finding using statistics of complexity and satellites |
| Synteny.nf | Alignment and synteny detection |
| nextflow.config | Enviromental excecution parameters for cluster |
| nextflow_ferropc.config | Enviromental excecution parameters for local cluster without queues and low resources |
| scripts/ | folder containing scripts for plotting and other small analysis |
| scripts/other.txt | text file containing small commands run on terminal for specific file conversions and opperations |
| example_sh | bash scripts to run the different nextflow pipelines, named acordingly | 
