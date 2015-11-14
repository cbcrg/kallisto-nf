Kallisto-NF
========

A nextflow implementation of Kallisto & Sleuth RNA-Seq Tools



Quick start 
-----------

Make sure you have all the required dependencies listed in the last section.

Install the Nextflow runtime by running the following command:

    $ curl -fsSL get.nextflow.io | bash


When done, you can launch the pipeline execution by entering the command shown below:

    $ nextflow run cbcrg/kallisto-nf
    

By default the pipeline is executed against the provided example dataset. 
Check the *Pipeline parameters*  section below to see how enter your data on the program 
command line.     
    


Pipeline parameters
-------------------

**--transcriptome**  
   
* The location of the transcriptome multi-fasta file. 
* It should end in '.fa' 
* Involved in the task: kallisto-index.
  * By default is set to the Kallisto-NF's localization: './tutorial/data/transcriptome/trascriptome.fa'
  `  $ nextflow run cbcrg/transcriptome-nf --transcriptome /home/user/my_transcriptome/example.fa  `
    
  
**--reads** 
   
* Specifies the location of the reads *fastq* file(s).
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)
* It must end in 'fastq'.
* Involved in the task: kallisto-mapping.
  * By default is set to the Kallisto-NF's location: './tutorial/data/*.fastq' 
  `  $ nextflow run cbcrg/kallisto-nf --reads '/home/dataset/*.fastq'`
* Read pairs can be specified using the a glob file patter like the following: 
  `  $ nextflow run cbcrg/kallisto-nf --reads '/home/dataset/*_{1,2}.fastq'`    
  

**--fragment_len**
* Specifies the average fragment length of the RNA-Seq library.
* This is required for mapping single-ended reads.
* Involved in the task: kallisto-mapping.
   * By default is set 180. 
  `  $ nextflow run cbcrg/kallisto-nf --fragment_len 180`


**--fragment_sd**
* Specifies the standard deviation of the fragment length in the RNA-Seq library.
* This is required for mapping single-ended reads.
* Involved in the task: kallisto-mapping.
   * By default this is set 20.  
  `  $ nextflow run cbcrg/kallisto-nf --fragment_sd 180`


**--bootstrap**
* Specifies the number of bootstrap samples for quantification of abundances.
* Involved in the task: kallisto-mapping.
   * By default this is set 100. 
  `  $ nextflow run cbcrg/kallisto-nf --bootstrap 100`


**--experiment**
* Specifies the location of the experimental design file. 
* The experimental design file provides Seulth with a link between the samples and the conditions for abundance testing. 
  * By default is set to the Kallisto-NF's location: './tutorial/experiment/high_seqinfo.txt'
  `  $ nextflow run cbcrg/kallisto-nf --experiment '/home/experiment/exp_design.txt'`


**--output** 
   
* Specifies the folder where the results will be stored for the user.  
* It does not matter if the folder does not exist.
  * By default is set to Kallisto-NF's folder: './results' 
  `  $ nextflow run cbcrg/kallisto-nf --output /home/user/my_results  `
  


Cluster support
---------------

Kallisto-NF execution relies on [Nextflow](http://www.nextflow.io) framework which provides an 
abstraction between the pipeline functional logic and the underlying processing system.

Thus it is possible to execute it on your computer or any cluster resource
manager without modifying it.

Currently the following clusters are supported:

  + Oracle/Univa/Open Grid Engine (SGE)
  + Platform LSF
  + SLURM
  + PBS/Torque


By default the pipeline is parallelized by spanning multiple threads in the machine where the script is launched.

To submit the execution to a SGE cluster create a file named `nextflow.config`, in the directory
where the pipeline is going to be launched, with the following content:

    process {
      executor='sge'
      queue='<your queue name>'
    }

In doing that, tasks will be executed through the `qsub` SGE command, and so your pipeline will behave like any
other SGE job script, with the benefit that *Nextflow* will automatically and transparently manage the tasks
synchronisation, file(s) staging/un-staging, etc.

Alternatively the same declaration can be defined in the file `$HOME/.nextflow/config`.

To lean more about the avaible settings and the configuration file read the Nextflow documentation 
 http://www.nextflow.io/docs/latest/config.html
  
  
Dependencies 
------------

 * Java 7+ 
 * Kallisto - http://pachterlab.github.io/kallisto/
 * Sleuth - http://pachterlab.github.io/sleuth/
 * R - https://www.r-project.org/
 
