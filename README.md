## *Crenicichla tuca* genome assembly
This branch contains nextflow scripts used for all assembly steps and analysis presented in *C. tuca* genome report.

<img src="https://github.com/mylena-s/assembly_project/blob/CtucaBranch/github.png?raw=true" width="1000">

**Reference**

All analysis were performed integrated in nextflow these workflows, except for:
* D-genies: run using an interactive local server.
* TOGA annotation: custom script runToga.sh
* FCS-NCBI: run on Usegalaxy.org
* OMark: run on https://omark.omabrowser.org/
* Small file manipulation scripts: python or bash scripts used for small taks can be found on /others

### Example scripts 
Examples of bash scripts used for executing the pipeline are found on /examples

All software was installed using mamba or otherwise stated in software/install.sh

# Scripts explanation
| Script | Uses |
-----------------
| readQC.nf | Read quality assessment before assembly |
| Assembly.nf | Genome assembly with ONT and HIFI reads |
| AssemblyCuration.nf | Scaffolding |
| AssemblyQC.nf | Assembly quality assessment, including BUSCO, kat, and other software |
| Annotation.nf | Gene and repeat annotation |	
| Centromere.nf | Centromere finding using statistics of complexity and satellites |
| Synteny.nf | Alignment and synteny detection |
| nextflow.config | Enviromental excecution parameters for cluster |
| nextflow_ferropc.config | Enviromental excecution parameters for local cluster without queues and low resources |
