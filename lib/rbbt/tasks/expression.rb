module CLSS

  input :cell_lines, :array, "Cell line names", %w(AGS SW620 DU145 Colo205 UACC62 SF295 A498 MDA-MB-468)
  task :progeny => :tsv do |cell_lines|
    ccle_cls = cell_lines.collect{|cl| ccle_cell_line(cl) }.compact
    tsv = TSV.open(RbbtMatrix.new(CCLE.gene_expression.find).to_average.data_file, :fields => ccle_cls)
    tsv.R <<-EOF
rbbt.require('progeny')
pathways = progeny(as.matrix(data), scale=FALSE)
data = pathways

    EOF
  end

end
