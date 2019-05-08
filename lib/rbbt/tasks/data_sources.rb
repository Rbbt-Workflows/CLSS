module CLSS
  helper :matrix_activity do |tsv_file,sample|
    m = RbbtMatrix.new tsv_file
    TmpFile.with_file do |outfile|
      tsv = begin
              TSV.open(m.to_name.to_average.to_activity.data_file, :type => :list, :cast => :to_f)
            rescue
              TSV.open(m.to_average.to_activity.data_file, :type => :list, :cast => :to_f)
            end
      tsv.each do |k,v|
        v_s = Misc.std_num_vector(v, -1, 1)
        v.replace(v_s.collect{|e| e.nan? ? nil : e })
      end
      raise ParameterException, "#{ sample } field not found in matrix" unless tsv.fields.include? sample
      tsv = tsv.slice(sample).transpose
      tsv.key_field = "id"
      tsv
    end
  end

  helper :achilles_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(Achilles.cell_line_info.tsv.keys).resolve(cell_line, nil, :max_candidates => 100, :threshold => -1000).first
  end

  helper :ccle_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(CCLE.cell_lines, :use_keys => true).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  helper :gdsc_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    Normalizer.new(GDSC.cell_lines, :use_keys => true).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  helper :mclp_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    cell_lines = MCLP.RPPA.tsv.fields
    Normalizer.new(cell_lines).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  helper :cosmic_cell_line do |cell_line|
    cl = GDSC.cell_lines.index(:target => "COSMIC cell line ID", :fields => ["Sample Name"], :persist => true)[cell_line]
    cl = GDSC.cell_line_equivalences.index(:target => "COSMIC cell line ID", :persist => true)[cell_line] if cl.nil?
    cl
  end

  helper :coread_cell_line do |cell_line|
    require 'rbbt/ner/rnorm'
    require 'rbbt/sources/COREAD_phospho_proteome'
    cell_lines = COREADPhosphoProteome.phosphosite_levels.tsv.fields
    Normalizer.new(cell_lines).resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
  end

  task :CCLE_to_GDSC => :tsv do
    require 'rbbt/ner/rnorm'
    norm = Normalizer.new(GDSC.cell_lines, :use_keys => true)
    tsv = TSV.setup({},"CCLE~GDSC#type=single")
    CCLE.cell_lines.tsv.keys.each do |cell_line|
      orig = cell_line
      cell_line = cell_line.split("_").first
      gdsc = norm.resolve(cell_line, nil, :max_candidates => 100, :threshold => 0).first
      tsv[orig] = gdsc
    end
    tsv
  end

  input :cell_line, :string, "Cell line name"
  input :cl_data_source, :select, "Source of cell line data", 'CCLE', :select_options => %w(CCLE GDSC)
  task :mRNA => :tsv do |cell_line,cl_data_source|
    case cl_data_source
    when "CCLE"
      ccle_cl = ccle_cell_line(cell_line)
      raise ParameterException, "CCLE cell line not found: " << cell_line if ccle_cl.nil?
      log :cell_line, "Using cell_line #{ccle_cl} (#{cell_line})"
      tsv = matrix_activity(CCLE.gene_expression.find, ccle_cl)
      current = tsv.keys.first
      tsv[cell_line] = tsv[current]
      tsv.delete(current) unless current == cell_line
    when "GDSC"
      gdsc_cl = gdsc_cell_line(cell_line)
      raise ParameterException, "GDSC cell line not found: " << cell_line if gdsc_cl.nil?
      log :cell_line, "Using cell_line #{gdsc_cl} (#{cell_line})"
      cosmic_cl = cosmic_cell_line(gdsc_cl)
      raise ParameterException, "COSMIC cell line not found: " << cell_line if cosmic_cl.nil?
      tsv = matrix_activity(GDSC.gene_expression.find, cosmic_cl)
      current = tsv.keys.first
      tsv[cell_line] = tsv[current]
      tsv.delete(current) unless current == cell_line
    end

    tsv.key_field = "id"
    tsv
  end
  
  input :cell_line, :string, "Cell line name"
  input :cl_data_source, :select, "Source of cell line data", 'CCLE', :select_options => %w(CCLE GDSC)
  task :genome => :tsv do |cell_line,cl_data_source|
    case cl_data_source
    when "CCLE"
      ccle_cl = ccle_cell_line(cell_line)
      raise ParameterException, "Cell line not found: " << cell_line if ccle_cl.nil?
      log :cell_line, "Using cell_line #{ccle_cl} (#{cell_line})"

      tsv = begin
              CCLE.gene_CNV.tsv(:type => :list, :fields => [ccle_cl]).change_key("Associated Gene Name", :identifiers  => Organism.identifiers(CCLE.organism)).transpose
            rescue
              raise ParameterException, "Cell line has no CNV data: " << cell_line 
            end

      current = tsv.keys.first

      tsv[cell_line] = tsv[current]
      tsv.delete(current) unless current == cell_line
    when "GDSC"
      gdsc_cl = gdsc_cell_line(cell_line)
      raise ParameterException, "Cell line not found: " << cell_line if gdsc_cl.nil?
      log :cell_line, "Using cell_line #{gdsc_cl} (#{cell_line})"

      tsv = begin
               tsv = GDSC.cell_line_CEFs.tsv
               cosmic_cl = cosmic_cell_line(gdsc_cl)
               raise ParameterException, "COSMIC cell line not found: " << cell_line if cosmic_cl.nil?
               gained = tsv[cosmic_cl].select{|c| c.include?("gain:") && c.include?("(") }.collect{|c| c.match(/\((.*)\)/)[1].split(",") }.flatten.collect{|g| g.strip}
               lost = tsv[cosmic_cl].select{|c| c.include?("loss:") && c.include?("(") }.collect{|c| c.match(/\((.*)\)/)[1].split(",") }.flatten.collect{|g| g.strip}
               tsv = TSV.setup({}, :key_field => "id", :fields => gained + lost, :type => :list)
               tsv[cell_line] = [1] * gained.length + [-1] * lost.length
               tsv
            rescue
               tsv = TSV.setup({}, :key_field => "id", :fields => [], :type => :list)
               tsv[cell_line] = []
               tsv
            end
    end
    tsv.key_field = "id"
    tsv
  end

  input :cell_line, :string, "Cell line name"
  input :binarize, :boolean, "Binarize RPPA values", true
  task :abundance => :tsv do |cell_line,binarize|
    mclp_cl = mclp_cell_line(cell_line)
    raise ParameterException, "Cell line not found: " << cell_line if mclp_cl.nil?

    log :cell_line, "Using cell_line #{mclp_cl} (#{cell_line})"
    if binarize
      tsv = matrix_activity(MCLP.RPPA.find, mclp_cl)
    else
      tsv = MCLP.RPPA.tsv
      tsv = tsv.slice(mclp_cl).transpose
      tsv.key_field = "id"
      tsv
    end

    annotations = Rbbt.share["MCLP_RPPA_annotation_20170222.tsv"].tsv :fields => ["Genes", "Correlation"], :type => :double

    abs = tsv.fields 
    all_genes = annotations.column("Genes").values.flatten.compact.reject{|g| g.empty? }.sort.uniq
    activity = TSV.setup({}, :key_field => "id", :fields => all_genes, :type => :double)
    tsv.each do |sample,values|
      res = [nil] * all_genes.length
      values.zip(abs).each do |value,ab|
        genes = annotations[ab]["Genes"]
        next if genes.empty?
        correlation = annotations[ab]["Correlation"].first.to_f
        next if correlation != 0
        genes.each do |gene|
          res[all_genes.index(gene)] = value
        end
      end
      activity.zip_new sample, res
    end

    activity[cell_line] = activity[mclp_cl]
    activity.delete mclp_cl unless mclp_cl == cell_line
    activity.key_field = "id"
    activity.to_list{|f| f.any? ? Misc.mean(f) : nil }
  end

  input :cell_line, :string, "Cell line name"
  input :binarize, :boolean, "Binarize RPPA values", true
  task :activity => :tsv do |cell_line,binarize|
    mclp_cl = mclp_cell_line(cell_line)
    raise ParameterException, "Cell line not found: " << cell_line if mclp_cl.nil?

    log :cell_line, "Using cell_line #{mclp_cl} (#{cell_line})"
    if binarize
      tsv = matrix_activity(MCLP.RPPA.find, mclp_cl)
    else
      tsv = MCLP.RPPA.tsv
      tsv = tsv.slice(mclp_cl).transpose
      tsv.key_field = "id"
      tsv
    end

    annotations = Rbbt.share["MCLP_RPPA_annotation_20170222.tsv"].tsv :fields => ["Genes", "Correlation"], :type => :double

    abs = tsv.fields 
    all_genes = annotations.column("Genes").values.flatten.compact.reject{|g| g.empty? }.sort.uniq
    activity = TSV.setup({}, :key_field => "id", :fields => all_genes, :type => :double)
    tsv.each do |sample,values|
      res = [nil] * all_genes.length
      values.zip(abs).each do |value,ab|
        gene = annotations[ab]["Genes"].first
        next if gene.nil? or gene.empty?
        correlation = annotations[ab]["Correlation"].first.to_f
        next if correlation == 0
        v = value.to_f * correlation
        res[all_genes.index(gene)] = v
      end
      activity.zip_new sample, res
    end

    activity[cell_line] = activity[mclp_cl]
    activity.delete mclp_cl unless mclp_cl == cell_line
    activity.key_field = "id"
    activity.to_list{|f| f.any? ? Misc.mean(f) : nil }
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
    raise parameterexception, "no cosmic cell line found in gdsc for cell line: " << cell_line if cosmic_cl.nil?
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
