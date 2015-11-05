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

import org.apache.commons.lang.StringUtils

/* 
 * Main Kallisto-NF pipeline script
 *
 * @authors
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */


params.transcriptome = "$baseDir/tutorial/transcriptome/Homo_sapiens.GRCh38.rel79.cdna.part.fa"
params.name          = "Expression Analysis"
params.primary       = "$baseDir/tutorial/reads/*.fastq"
params.secondary     = null
params.exp           = "$baseDir/tutorial/experiment/hiseq_info.txt"
params.db            = "db/"
params.mapper        = "kallisto"
params.cpus          = 1
params.output        = "results/"


log.info "K A L L I S T O - N F  ~  version 0.1"
log.info "================================="
log.info "name                   : ${params.name}"
log.info "transcriptome          : ${params.transcriptome}"
log.info "primary                : ${params.primary}"
log.info "secondary              : ${params.secondary}"
log.info "experimental design    : ${params.exp}"
log.info "database path          : ${params.db}"
log.info "output                 : ${params.output}"
log.info "mapper                 : ${params.mapper}"
log.info "\n"


/*
 * Input parameters validation
 */

if( !(params.mapper in ['kallisto'])) { exit 1, "Invalid mapper tool: '${params.mapper}'" }

transcriptome_file     = file(params.transcriptome)
primary_reads          = files(params.primary).sort()

if (params.secondary != null) {
    secondary_reads  = files(params.secondary).sort()
}
else {
    secondary_reads = []
}  

exp_file               = file(params.exp) 
result_path            = file(params.output)
dbPath                 = file(params.db)

/*
 * validate input files
 */
if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"

if( !exp_file.exists() ) exit 1, "Missing experimental design file: ${exp_file}"

if( !result_path.exists() && !result_path.mkdirs() ) {
    exit 3, "Cannot create output folder: $result_path -- Check file system access permission"
}

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
    
    storeDir { dbPath }

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
    scratch false

    storeDir { dbPath }
    
    input:
    file transcriptome_index from transcriptome_index.first()
    file primary_reads
    val read_names

    output:
    file "kallisto/${read_names}/abundance.h5" into kallisto_h5
    file "kallisto/${read_names}/abundance.tsv" into kallisto_tsv
    file "kallisto/${read_names}/run_info.json" into kallisto_run_info   

    script:
    //
    // Kallisto tools mapper
    //

    if( params.mapper == 'kallisto' &&  secondary_reads.size() == 0) {
        """
        test -d "kallisto/" || mkdir -p "kallisto/"
 
        test -d "${result_path}/kallisto/" || mkdir -p "${result_path}/kallisto/"
 
        kallisto quant --single -l 200 -s 3 -b 100 -t 8 -i ${transcriptome_index} -o kallisto/${read_names} ${primary_reads} 

        cp -r kallisto/${read_names} "${result_path}/kallisto/." 
        """
    }  
    else {
        """
        kallisto quant -i transcriptome.index -o ${read_names} ${primary_reads} ${secondary_reads}
        """
    }

}


process sleuth {
 
    executor='local'   

    input:
    file 'h5_*' from kallisto_h5.toList()
    file 'tsv_*' from kallisto_tsv.toList()
    file 'run_info_*' from kallisto_run_info.toList()


    script:
    //
    // Setup sleuth R dependancies and environment
    //
 
    """
    #!/usr/bin/env /nfs/users/cn/efloden/R-3.2.2/bin/Rscript
    library("sleuth")

    sample_id <- dir(file.path("$baseDir","results", "kallisto"))
    sample_id

    kal_dirs <- sapply(sample_id, function(id) file.path("$baseDir", "results", "kallisto", id))

    s2c <- read.table(file.path("$baseDir","data", "experiments", "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
    s2c <- dplyr::select(s2c, sample = run_accession, condition)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)

    mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

    so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g)
    so <- sleuth_fit(so)
    so <- sleuth_wt(so, 'conditionscramble')

    sleuth_live(so)
    """
}

// ===================== UTILITY FUNCTIONS ============================


/*
 * Given a path returns a sorted list files matching it.
 * The path can contains wildcards characters '*' and '?'
 */
def List<Path> findReads( String fileName ) {
    def result = []
    if( fileName.contains('*') || fileName.contains('?') ) {
        def path = file(fileName)
        def parent = path.parent
        def filePattern = path.getName().replace("?", ".?").replace("*", ".*")
        parent.eachFileMatch(~/$filePattern/) { result << it }
        result = result.sort()
    }
    else {
        result << file(fileName)
    }

    return result
}

def bestMatch( Path file1, Path file2, boolean singlePair = true, boolean pairedEnd = false) {
    bestMatch( file1.baseName, file2.baseName, singlePair, pairedEnd)
}

def bestMatch( String n1, String n2, boolean singlePair = true, boolean pairedEnd = false) {

    def index = StringUtils.indexOfDifference(n1, n2)

    if (!pairedEnd) {
        String match = n1
        match = trimReadName(match)
        return [match, null]
    }

  
    else {
        if( !singlePair ) {
            if( index == -1 ) {
                // this mean the two file names are identical, something is wrong
                return [null, "Missing entry for read pair: '$n1'"]
            }
            else if( index == 0 ) {
                // this mean the two file names are completely different
                return [null, "Not a valid read pair -- primary: $n1 ~ secondary: $n2"]
            }
        }

        String match = index ? n1.subSequence(0,index) : n1
        match = trimReadName(match)
        if( !match ) {
            return [null, "Missing common name for read pair -- primary: $n1 ~ secondary: $n2 "]
        }

        return [match, null]
    }
}

def trimReadName( Path file1 ) {
    trimReadName( file1.baseName ) 
}
def trimReadName( String name ) {
    name.replaceAll(/^[^a-zA-Z0-9]*/,'').replaceAll(/[^a-zA-Z0-9]*$/,'')
}


// ===================== UNIT TESTS ============================

def testFindReads() {

    def path = File.createTempDir().toPath()
    try {
        def file1 = path.resolve('alpha_1.fastq'); file1.text = 'file1'
        def file2 = path.resolve('alpha_2.fastq'); file2.text = 'file2'
        def file3 = path.resolve('gamma_1.fastq'); file3.text = 'file3'
        def file4 = path.resolve('gamma_2.fastq'); file4.text = 'file4'

        assert files("$path/alpha_1.fastq") == [file1]
        assert files("$path/*_1.fastq") == [file1, file3]
        assert files("$path/*_2.fastq") == [file2, file4]
    }
    finally {
        path.deleteDir()
    }

}

def testTrimReadName() {
    assert trimReadName('abc') == 'abc'
    assert trimReadName('a_b_c__') == 'a_b_c'
    assert trimReadName('__a_b_c__') == 'a_b_c'
}

def testBestMach() {

    assert bestMatch('abc_1', 'abc_2') == ['abc', null]
    assert bestMatch('aaa', 'bbb') == ['aaa', null]
    assert bestMatch('_', 'bbb') == [null, "Missing common name for read pair -- primary: _ ~ secondary: bbb "]

    assert bestMatch('abc_1', 'abc_2', false) == ['abc', null]
    assert bestMatch('aaa', 'bbb', false) == [null, 'Not a valid read pair -- primary: aaa ~ secondary: bbb' ]

}
