shell.executable("/bin/bash")

import pandas
import altair

IDS, = glob_wildcards("genomes/{id}.fa")

rule all:
	input: expand("hists/{sample}-ratioplot-full.html", sample=IDS)

rule prodigal:
	input: "genomes/{id}.fa"
	output: "proteins/{id}.faa"
	shell: "prodigal -a {output} -q -i {input}"

rule diamond:
	input: "proteins/{id}.faa"
	output: "lengths/{id}.data"
	threads: 4
	params:
		db="uniprot_trembl.diamond.dmnd",
		of="6 qlen slen"
	shell: "diamond blastp --threads {threads} --max-target-seqs 1 --db {params.db} --query {input} --outfmt {params.of} --out {output}"

rule hist:
	input: "lengths/{id}.data"
	output: "hists/{id}-ratioplot-full.html"
	#	shell: "scripts/hist.R {input} {output}"
	run:
		df = pandas.read_csv(str(input), sep = "\t", names = ["qlen", "slen"])
	    # the 'slen' value is always one less than the query because database doesn't include stop codons
		df['slen'] = df['slen'] + 1

	    # add a value for the query length divided by the sequence length)
		df['codingRatio'] = df['qlen'] / df['slen']

		# add a 'gene number', assumes proteins were searched in order
		df['geneNumber'] = df.index + 1

		hist = altair.Chart(df)\
    	.mark_bar(clip = True)\
        .encode(x = altair.X('codingRatio',
                             bin = altair.Bin(step = 0.1),
                             scale = altair.Scale(domain=(0, 2)),
                             axis = altair.Axis(title='Query/Reference Ratio')
                            ),
                y = 'count()',
               tooltip = 'count()')\
        .configure_mark(
            fill = 'red',
            stroke = 'black')

		histzoom = altair.Chart(df)\
		.mark_bar(clip = True)\
		.encode(x = altair.X('codingRatio',
	                         bin = altair.Bin(step = 0.1),
	                         scale = altair.Scale(domain=(0, 2)),
	                         axis = altair.Axis(title='Query/Reference Ratio')
	                        ),
	            y = altair.Y('count()',
	                        scale = altair.Scale(domain = (0,100))
	                        ),
	           tooltip = 'count()')\
	    .configure_mark(
	        fill = ' #c658dd ',
	        stroke = 'black')

		genomeRatio = altair.Chart(df)\
	    .mark_line(clip = True)\
	    .encode(x = altair.X('geneNumber',
	                         scale = altair.Scale(domain = (0, len(df.index)))
	                        ),
	            y = altair.Y('codingRatio',
	                        scale = altair.Scale(type = 'log')),
	            tooltip = 'codingRatio')\
	    .interactive()

	    # save outputs
		#hist.save(str(output))
		output = str(output)
		hist.save(output)
		histzoom.save(output.replace('-ratioplot-full.html', '-ratioplot-zoom.html'))
		genomeRatio.save(output.replace('-ratioplot-full.html', '-ratioplot-genome.html'))
