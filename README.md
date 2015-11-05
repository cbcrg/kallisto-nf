Kallisto-NF
========

A nextflow implementation of Kallisto & Slueth RNA-Seq Tools



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
   
* The location of the genome multi-fasta file. 
* It should end in '.fa' 
* Involved in the task: index.
  * By default is set to the Kallisto-NF's localization: './tutorial/data/transcriptome/trascriptome.fa'
  `  $ nextflow run cbcrg/transcriptome-nf --transcriptome /home/user/my_fastas/example.fa  `
    
  
**--primary** 
   
* Specifies the location of the primary reads *fastq* file.
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)
* It must end in 'fastq'.
* Involved in the task: rna-pipeline.
  * By default is set to the Kallisto-NF's location: './tutorial/data/test_1.fastq' 
  `  $ nextflow run cbcrg/kallisto-nf --primary '/home/dataset/*_1.fastq'`
  
  
**--secondary** 
   
* Specifies the location of the secondary reads *fastq* file if paired end data are used.
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
   value by single quote characters (see the example below)
* It must end in '_2.fastq'.  
* Involved in the task: rna-pipeline.  
  * By default is set to the Kallisto-NF's location: './tutorial/data/test_2.fastq' 
  `  $ nextflow run cbcrg/kallisto-nf --secondary '/home/dataset/*_2.fastq'`

**--exp**
* Specifies the location of the experimental design file

**--cpus** 
   
* Sets the number of CPUs used in every tasks (default 1).  
* Involved in the task of kallisto bootstraps
  * By default is set to the number of the available cores.  
  `  $ nextflow run cbcrg/kallisto-nf --cpus 10  `
  
  
**--output** 
   
* Specifies the folder where the results will be stored for the user.  
* It does not matter if the folder does not exist.
  * By default is set to Grape-NF's folder: './results' 
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
 * Sleuth- http://pachterlab.github.io/sleuth/
 * R - https://www.r-project.org/
 
