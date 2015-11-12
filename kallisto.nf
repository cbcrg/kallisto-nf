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
params.primary       = "$baseDir/tutorial/reads/*.fastq"
params.secondary     = null
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '100'
params.experiment    = "$baseDir/tutorial/experiment/hiseq_info.txt"
params.mapper        = "kallisto"
params.output        = "results/"


log.info "K A L L I S T O - N F  ~  version 0.2"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "transcriptome          : ${params.transcriptome}"
log.info "primary                : ${params.primary}"
log.info "secondary              : ${params.secondary}"
log.info "fragment length        : ${params.fragment_len} nt"
log.info "fragment SD            : ${params.fragment_sd} nt"
log.info "bootstraps             : ${params.bootstrap}"
log.info "experimental design    : ${params.experiment}"
log.info "output                 : ${params.output}"
log.info "mapper                 : ${params.mapper}"
log.info "\n"


/*
 * Input parameters validation
 */

if( !(params.mapper in ['kallisto'])) { exit 1, "Invalid mapper tool: '${params.mapper}'" }

transcriptome_file     = file(params.transcriptome)
primary_reads          = files(params.primary).sort()
secondary_reads = params.secondary ?  files(params.secondary).sort() :  []
exp_file               = file(params.experiment) 
result_path            = file(params.output)

/*
 * validate input files
 */
if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"

if( !exp_file.exists() ) exit 1, "Missing experimental design file: ${exp_file}"

read_names = []
for( file in primary_reads ) {
    String name = trimReadName(file)
    read_names << name
}
read_names.sort()


log.info "Read pairs: $read_names"
log.debug "Primary reads: ${primary_reads *. name }"
if (secondary_reads.size() != 0) {
  log.debug "Secondary reads: ${secondary_reads *. name }"
}

process index {

    input:
    file transcriptome_file
    
    output:
    file "transcriptome.index" into transcriptome_index
      
    script:
    //
    // Kallisto tools mapper index
    //
    if( params.mapper=='kallisto' )
        """
        kallisto index -i transcriptome.index ${transcriptome_file}
        """

}


process mapping {
    publishDir result_path

    input:
    file transcriptome_index from transcriptome_index.first()
    file primary_reads
    val read_names

    output:
    file "kallisto_${read_names}" into kallisto_out  

    script:
    //
    // Kallisto tools mapper
    //

    if( params.mapper == 'kallisto' &&  secondary_reads.size() == 0) {
        """
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} -t 8 -i ${transcriptome_index} -o kallisto_${read_names} ${primary_reads} 
        """
    }  
    else {
        """
        kallisto quant -i transcriptome.index -o ${read_names} ${primary_reads} ${secondary_reads}
        """
    }

}


process sleuth {
    publishDir result_path
 
    input:
    file all_files from kallisto_out.toList()
    file 'hiseq_info.txt' from exp_file

    output: 
    file 'sleuth_object.so'

    script:
    //
    // Setup sleuth R dependancies and environment
    //
 
    """
    mkdir kallisto
    mv $all_files kallisto 
    sleuth.R kallisto hiseq_info.txt
    """
}


// ===================== UTILITY FUNCTIONS ============================


def trimReadName( Path file1 ) {
    trimReadName( file1.baseName ) 
}
def trimReadName( String name ) {
    name.replaceAll(/^[^a-zA-Z0-9]*/,'').replaceAll(/[^a-zA-Z0-9]*$/,'')
}


// ===================== UNIT TESTS ============================


def testTrimReadName() {
    assert trimReadName('abc') == 'abc'
    assert trimReadName('a_b_c__') == 'a_b_c'
    assert trimReadName('__a_b_c__') == 'a_b_c'
}


