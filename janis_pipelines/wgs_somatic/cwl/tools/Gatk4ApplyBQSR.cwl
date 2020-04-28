#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: Gatk4ApplyBQSR
doc: |-
  Apply base quality score recalibration: This tool performs the second pass in a two-stage 
  process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the 
  base qualities of the input reads based on the recalibration table produced by the 
  BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.

  Summary of the BQSR procedure: The goal of this procedure is to correct for systematic bias 
  that affect the assignment of base quality scores by the sequencer. The first pass consists 
  of calculating error empirically and finding patterns in how error varies with basecall 
  features over all bases. The relevant observations are written to a recalibration table. 
  The second pass consists of applying numerical corrections to each individual basecall 
  based on the patterns identified in the first step (recorded in the recalibration table) 
  and write out the recalibrated data to a new BAM or CRAM file.

  - This tool replaces the use of PrintReads for the application of base quality score 
      recalibration as practiced in earlier versions of GATK (2.x and 3.x).
  - You should only run ApplyBQSR with the covariates table created from the input BAM or CRAM file(s).
  - Original qualities can be retained in the output file under the "OQ" tag if desired. 
      See the `--emit-original-quals` argument for details.
requirements:
  DockerRequirement:
    dockerPull: broadinstitute/gatk:4.1.3.0
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
inputs:
- id: bam
  label: bam
  doc: The SAM/BAM/CRAM file containing reads.
  type: File
  secondaryFiles: |-
    ${

            function resolveSecondary(base, secPattern) {
              if (secPattern[0] == "^") {
                var spl = base.split(".");
                var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
              }
              return base + secPattern
            }

            return [
                    {
                        location: resolveSecondary(self.location, "^.bai"),
                        basename: resolveSecondary(self.basename, ".bai"),
                        class: "File",
                    }
            ];

    }
  inputBinding:
    prefix: -I
    position: 10
- id: reference
  label: reference
  doc: Reference sequence
  type: File
  secondaryFiles:
  - .fai
  - .amb
  - .ann
  - .bwt
  - .pac
  - .sa
  - ^.dict
  inputBinding:
    prefix: -R
- id: outputFilename
  label: outputFilename
  doc: Write output to this file
  type: string
  default: generated.bam
  inputBinding:
    prefix: -O
- id: recalFile
  label: recalFile
  doc: Input recalibration table for BQSR
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --bqsr-recal-file
- id: intervals
  label: intervals
  doc: -L (BASE) One or more genomic intervals over which to operate
  type:
  - File
  - 'null'
  inputBinding:
    prefix: --intervals
- id: tmpDir
  label: tmpDir
  doc: Temp directory to use.
  type: string
  default: /tmp/
  inputBinding:
    prefix: --tmp-dir
    position: 11
outputs:
- id: out
  label: out
  type: File
  secondaryFiles: |-
    ${

            function resolveSecondary(base, secPattern) {
              if (secPattern[0] == "^") {
                var spl = base.split(".");
                var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
              }
              return base + secPattern
            }
            return [
                    {
                        path: resolveSecondary(self.path, "^.bai"),
                        basename: resolveSecondary(self.basename, ".bai"),
                        class: "File",
                    }
            ];

    }
  outputBinding:
    glob: $(inputs.outputFilename)
baseCommand:
- gatk
- ApplyBQSR
id: Gatk4ApplyBQSR
