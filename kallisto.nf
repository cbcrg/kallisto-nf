/*
 * Copyright (c) 2015, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'Kallisto-NF'.
 *
 *   Kallisto-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Kallisto-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Kallisto-NF.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main Kallisto-NF pipeline script
 *
 * @authors
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */


params.transcriptome = "$baseDir/tutorial/transcriptome/transcriptome.fa"
params.name          = "RNA-Seq Abundance Analysis"
params.reads         = "$baseDir/tutorial/reads/*.fastq"
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '100'
params.threads       = '1'
params.experiment    = "$baseDir/tutorial/experiment/hiseq_info.txt"
params.output        = "results/"


log.info "K A L L I S T O - N F  ~  version 0.3"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "reads                  : ${params.reads}"
log.info "transcriptome          : ${params.transcriptome}"
log.info "fragment length        : ${params.fragment_len} nt"
log.info "fragment SD            : ${params.fragment_sd} nt"
log.info "bootstraps             : ${params.bootstrap}"
log.info "threads                : ${params.threads}"
log.info "experimental design    : ${params.experiment}"
log.info "output                 : ${params.output}"
log.info "\n"


/*
 * Input parameters validation
 */

transcriptome_file     = file(params.transcriptome)
exp_file               = file(params.experiment) 

/*
 * validate input files
 */
if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"

if( !exp_file.exists() ) exit 1, "Missing experimental design file: ${exp_file}"

/*
 * Create a channel for read files 
 */
 
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path -> 
       def prefix = readPrefix(path, params.reads)
       tuple(prefix, path) 
    }
    .groupTuple(sort: true)
    .set { read_files } 


process index {
    input:
    file transcriptome_file
    
    output:
    file "transcriptome.index" into transcriptome_index
      
    script:
    //
    // Kallisto tools mapper index
    //
    """
    kallisto index -i transcriptome.index ${transcriptome_file}
    """
}


process mapping {
    tag "reads: $name"

    input:
    file transcriptome_index from transcriptome_index.first()
    set val(name), file(reads:'*') from read_files

    output:
    file "kallisto_${name}" into kallisto_out_dirs 

    script:
    //
    // Kallisto tools mapper
    //
    def single = reads instanceof Path
    if( !single ) {
        """
        mkdir kallisto_${name}
        kallisto quant -b ${params.bootstrap} -i transcriptome.index -t ${params.threads} -o kallisto_${name} ${reads}
        """
    }  
    else {
        """
        mkdir kallisto_${name}
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} -i ${transcriptome_index} -t ${params.threads} -o kallisto_${name} ${reads}
        """
    }

}


process sleuth {
    input:
    file 'kallisto/*' from kallisto_out_dirs.toSortedList()   
    file exp_file

    output: 
    file 'sleuth_object.so'
    file 'gene_table_results.txt'

    script:
    //
    // Setup sleuth R dependancies and environment
    //
 
    """
    sleuth.R kallisto ${exp_file}
    """
}


// ===================== UTILITY FUNCTIONS ============================


/* 
 * Helper function, given a file Path 
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 * 
 * For example: 
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 * 
 * Returns: 
 *   'file_alpha'
 */
 
def readPrefix( Path actual, template ) {

    final fileName = actual.getFileName().toString()

    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') ) 
        filePattern = '*' + filePattern 
  
    def regex = filePattern
                    .replace('.','\\.')
                    .replace('*','(.*)')
                    .replace('?','(.?)')
                    .replace('{','(?:')
                    .replace('}',')')
                    .replace(',','|')

    def matcher = (fileName =~ /$regex/)
    if( matcher.matches() ) {  
        def end = matcher.end(matcher.groupCount() )      
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') ) 
          prefix=prefix[0..-2]
          
        return prefix
    }
    
    return fileName
}


