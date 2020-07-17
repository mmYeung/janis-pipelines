version development

task bedtoolsgenomeCoverageBed {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    Int? runtime_seconds
    Int? runtime_disks
    Boolean? depth
    Boolean? depthZero
    Boolean? BedGraphFormat
    Boolean? BedGraphFormata
    Boolean? split
    String? strand
    Boolean? pairEnd
    Boolean? fragmentSize
    Boolean? du
    Boolean? fivePos
    Boolean? threePos
    Int? max
    Float? scale
    Boolean? trackline
    String? trackopts
    File? inputBam
    File? inputBed
    File? inputFile
    File? genome
  }
  command <<<
    genomeCoverageBed \
      ~{if defined(depth) then "-d" else ""} \
      ~{if defined(depthZero) then "-dz" else ""} \
      ~{if defined(BedGraphFormat) then "-bg" else ""} \
      ~{if defined(BedGraphFormata) then "-bga" else ""} \
      ~{if defined(split) then "-split" else ""} \
      ~{if defined(strand) then ("-strand '" + strand + "'") else ""} \
      ~{if defined(pairEnd) then "-pc" else ""} \
      ~{if defined(fragmentSize) then "-fs" else ""} \
      ~{if defined(du) then "-du" else ""} \
      ~{if defined(fivePos) then "-5" else ""} \
      ~{if defined(threePos) then "-3" else ""} \
      ~{if defined(max) then ("-max " + max) else ''} \
      ~{if defined(scale) then ("-scale " + scale) else ''} \
      ~{if defined(trackline) then "-trackline" else ""} \
      ~{if defined(trackopts) then ("-trackopts '" + trackopts + "'") else ""} \
      ~{if defined(inputBam) then ("-ibam '" + inputBam + "'") else ""} \
      ~{if defined(inputBed) then ("-iBed '" + inputBed + "'") else ""} \
      ~{if defined(inputFile) then ("-i '" + inputFile + "'") else ""} \
      ~{if defined(genome) then ("-g '" + genome + "'") else ""}
  >>>
  runtime {
    cpu: select_first([runtime_cpu, 1])
    disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
    docker: "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    duration: select_first([runtime_seconds, 86400])
    memory: "~{select_first([runtime_memory, 8, 4])}G"
    preemptible: 2
  }
  output {
    File out = stdout()
  }
}