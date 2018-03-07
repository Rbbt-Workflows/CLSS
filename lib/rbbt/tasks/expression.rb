module CLSS

  input :cell_lines, :array, "Cell line names", %w(AGS SW620 DU145 Colo205 UACC62 SF295 A498 MDA-MB-468)
  input :dataset, :select, "Expression dataset", 'CCLE', :select_options => %w(CCLE GDSC)
  task :progeny => :tsv do |cell_lines,dataset|
    tsv = nil
    if dataset == "CCLE"
      ccle_cls = cell_lines.collect{|cl| ccle_cell_line(cl) }.compact
      tsv = TSV.open(RbbtMatrix.new(CCLE.gene_expression.find).to_average.data_file, :fields => ccle_cls)
    else
      gdsc_cls = cell_lines.collect{|cl| gdsc_cell_line(cl) }.compact
      name2cosmic = GDSC.cell_lines.index :target => "COSMIC cell line ID"
      cosmic2name = GDSC.cell_lines.index :target => "Sample Name"
      gdsc_cls = name2cosmic.values_at *gdsc_cls
      tsv = TSV.open(RbbtMatrix.new(GDSC.gene_expression.find).to_gene.data_file, :fields => gdsc_cls).change_key("Associated Gene Name", :identifiers => Organism.identifiers(GDSC.organism))
      tsv.key_field = "COSMIC cell line ID"
      tsv.fields = cosmic2name.values_at *tsv.fields
    end

    tsv.R <<-EOF
rbbt.require('progeny')
pathways = progeny(as.matrix(data), scale=TRUE)
data = pathways
    EOF
  end

end

