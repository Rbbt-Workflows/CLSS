require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/matrix'
require 'rbbt/matrix/barcode'

require 'rbbt/sources/MCLP'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/CLSS'

Workflow.require_workflow "CCLE"
Workflow.require_workflow "GDSC"
Workflow.require_workflow "Achilles"

module CLSS
  extend Workflow

  helper :matrix_activity do |tsv_file,sample|
    m = RbbtMatrix.new tsv_file
    TmpFile.with_file do |outfile|
      tsv = TSV.open(m.to_average.to_activity.data_file, :type => :list, :cast => :to_f)
      tsv.each do |k,v|
        v_s = Misc.std_num_vector(v, -1, 1)
        v.replace(v_s.collect{|e| e.nan? ? nil : e })
      end
      tsv = tsv.transpose.select(sample)
      tsv.key_field = "id"
      tsv
    end
  end

  helper :achilles_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(Achilles.cell_line_info.tsv.keys).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  helper :ccle_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(CCLE.cell_lines).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  helper :gdsc_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(GDSC.cell_lines).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  helper :mclp_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    cell_lines = MCLP.RPPA.tsv.fields
    Normalizer.new(cell_lines).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  helper :cosmic_cell_line do |cell_line|
    GDSC.cell_line_equivalences.index(:target => "COSMIC cell line ID")[cell_line]
  end

  input :cell_line, :string, "Cell line name"
  task :mRNA => :tsv do |cell_line|
    ccle_cl = ccle_cell_line(cell_line)
    raise ParameterException, "Cell line not found: " << cell_line if ccle_cl.nil?
    log :cell_line, "Using cell_line #{ccle_cl} (#{cell_line})"
    tsv = matrix_activity(CCLE.gene_expression.find, ccle_cl)
    current = tsv.keys.first
    tsv[cell_line] = tsv[current]
    tsv.delete(current)
    tsv.key_field = "id"
    tsv
  end
  
  input :cell_line, :string, "Cell line name"
  task :genome => :tsv do |cell_line|
    ccle_cl = ccle_cell_line(cell_line)
    raise ParameterException, "Cell line not found: " << cell_line if ccle_cl.nil?
    log :cell_line, "Using cell_line #{ccle_cl} (#{cell_line})"

    tsv = CCLE.gene_CNV.tsv(:type => :list, :fields => [ccle_cl]).change_key("Associated Gene Name", :identifiers  => Organism.identifiers(CCLE.organism)).transpose
    current = tsv.keys.first
    tsv[cell_line] = tsv[current]
    tsv.delete(current)
    tsv.key_field = "id"
    tsv
  end

  input :cell_line, :string, "Cell line name"
  task :activity => :tsv do |cell_line|
    mclp_cl = mclp_cell_line(cell_line)
    raise ParameterException, "Cell line not found: " << cell_line if mclp_cl.nil?

    log :cell_line, "Using cell_line #{mclp_cl} (#{cell_line})"
    tsv = matrix_activity(MCLP.RPPA.find, mclp_cl)

    annotations = Rbbt.share["MCLP_RPPA_annotation_20170222.tsv"].tsv :fields => ["Genes", "Correlation"], :type => :double

    abs = tsv.fields 
    genes = annotations.column("Genes").values.flatten.compact.reject{|g| g.empty? }.sort.uniq
    activity = TSV.setup({}, :key_field => "id", :fields => genes)
    tsv.each do |sample,values|
      res = [nil] * genes.length
      values.zip(abs).each do |value,ab|
        gene = annotations[ab]["Genes"].first
        next if gene.nil? or gene.empty?
        correlation = annotations[ab]["Correlation"].first.to_f
        v = value.to_f * correlation
        res[genes.index(gene)] = v
      end
      activity[sample] = res
    end

    activity[cell_line] = activity[mclp_cl]
    activity.delete mclp_cl unless mclp_cl == cell_line
    activity
  end

  input :cell_line, :string, "Cell line name"
  input :ceres_threshold, :float, "Threshold for selection of essential genes from Achilles", -0.8
  task :achilles_essential_genes => :array do |cell_line,thr|
    acl = achilles_cell_line(cell_line)
    raise ParameterException, "Cell line not found: " << cell_line if acl.nil?

    data = Achilles.ceres_gene_effects.tsv :fields => [acl], :type => :single, :cast => :to_f
    data.select{|k,value| value < thr }.keys
  end


  dep :activity
  task :steady_states_activity => :tsv do
    tsv = step(:activity).load.to_list.transpose.to_single
    tsv.process tsv.fields.first do |v|
      v = v.to_f
      if v.nil?
        "-"
      elsif v > 0.3
        "1"
      elsif v < -0.3
        "0"
      else
        "-"
      end
    end
  end

  input :cell_line, :string, "Cell line name"
  task :ccle_mutations => :array do |cell_line|
    ccle_cl = ccle_cell_line(cell_line)
    tsv = CCLE.oncomap_maf.tsv(:fields => ["Associated Gene Name"], :type => :flat)
    tsv[ccle_cl]
  end

  input :cell_line, :string, "Cell line name"
  task :gdsc_mutations => :array do |cell_line|
    tsv = GDSC.cell_line_CEFs.tsv
    cosmic_cl = cosmic_cell_line(cell_line)
    raise ParameterException, "No COSMIC cell line found in GDSC for cell line: " << cell_line if cosmic_cl.nil?
    tsv[cosmic_cl].select{|c| c.include? "_mut" }.collect{|c| c.split("_").first}.uniq
  end

  dep :gdsc_mutations
  task :drug_mutations => :tsv do
    stream = StringIO.new
    stream << "#Drug\tTarget" << "\n"
    TSV.traverse Rbbt.share["drugs.tab"], :type => :flat, :into => stream do |d,v|
      [("Drug_" << d.first.gsub(/\s+/,'_')),["inhibits"] + v].flatten * "\t" + "\n"
    end

    mutations = step(:gdsc_mutations).load


    #oncogenes = Rbbt.share.oncogenes.list
    
    Workflow.require_workflow "Genomics"
    require 'rbbt/sources/oncodrive_role'
    require 'rbbt/entity/gene'
    oncogenes = OncoDriveROLE.oncogenes.name


    inhibited = []
    activated = []
    TSV.traverse mutations, :type => :array do |mutation|
      gene = mutation.split("_").first
      if oncogenes.include?(gene)
        activated << gene
      else
        inhibited << gene
      end
    end

    stream << ["Mutation_activations", ["activates", activated]].flatten * "\t" << "\n"
    stream << ["Mutation_inhibitions", ["inhibits", inhibited]].flatten * "\t" << "\n"

    stream.rewind
    stream.read
  end


  dep :drug_mutations
  task :perturbations => :tsv do

    drugs = []
    mutations = []
    TSV.traverse step(:drug_mutations), :type => :array do |line|
      next if line =~ /^#/
      alt = line.split("\t").first
      if alt =~ /Drug/
        drugs << alt
      else
        mutations << alt
      end
    end
    str = ""

    drugs.each do |drug|
      str << drug << "\n"
    end
    mutations.each do |m|
      str << m << "\n"
    end

    str << [mutations] * "\t" << "\n"

    drugs.each do |drug|
      str << (mutations + [drug]) * "\t" << "\n"
      str << ([mutations.first] + [drug]) * "\t" << "\n"
      str << ([mutations.last] + [drug]) * "\t" << "\n"
    end

    str
  end
end

require 'rbbt/tasks/paradigm'

#require 'CLSS/tasks/basic.rb'

#require 'rbbt/knowledge_base/CLSS'
#require 'rbbt/entity/CLSS'

