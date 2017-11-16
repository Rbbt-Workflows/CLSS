Workflow.require_workflow "Paradigm"
Workflow.require_workflow "FNL"
Workflow.require_workflow "Viper"
Workflow.require_workflow "ROMA"

module CLSS

  dep FNL, :regulon, :jobname => "Default"
  input :st_pathway, :text, "Signal transduction pathway in Paradigm", nil, :required => true, :stream => false
  task :st_tf_pathway => :text do |pth|
    pth = pth.read if IO === pth or File === pth
    pth = Misc.is_filename?(pth)? Open.read(pth) : pth
    lines = pth.split("\n")
    proteins = lines.collect{|l| l =~ /^protein\t(.*)/ and $1 }.compact

    require 'rbbt/sources/organism'
    require 'rbbt/sources/signor'
    signor = Signor.protein_protein.tsv
    signor.identifiers = Organism.identifiers(Signor.organism)
    signor.key_field = "UniProt/SwissProt Accession"
    signor = signor.change_key "Associated Gene Name"
    signor = signor.swap_id "Target (UniProt/SwissProt Accession)", "Associated Gene Name"
    i = Association.index(signor, :persist_file => file('signor.index'))

    regulon = step(:regulon).load 
    tfs = regulon.keys

    ms = i.subset(proteins, tfs)
    i.unnamed = false
    new_prots = Set.new
    new_assocs = []
    ms.each do |_m|
      info = i[_m]
      tf, tg = _m.split("~")
      effect = info["Effect"]

      next if effect.include? 'transcription'

      sign = nil
      sign = 1 if effect.split(/[^\w]/).collect{|w| w.downcase}.include? 'up'
      sign = -1 if effect.split(/[^\w]/).collect{|w| w.downcase}.include? 'down'
      next if sign.nil?

      new_prots << tf
      new_prots << tg

      new_assocs << [tf, tg, sign]
    end

    new_lines = []
    (new_prots - proteins).each do |p|
      new_lines << "protein\t#{p}"
    end
    new_lines.concat lines

    new_assocs.each do |a|
      tf, tg, sign = a
      symbol = sign == 1 ? '-a>' : '-a|'
      new_lines << [tf, tg, symbol] * "\t"
    end

    new_lines * "\n" + "\n"
  end

  dep FNL, :regulon, :compute => :produce, :jobname => "Default"
  task :regulon_modules => :text do
    TSV.traverse step(:regulon), :into => :stream do |tf, values|
      tgs = Misc.zip_fields(values).collect{|tf,s|s = 1 if (s.nil? or s.empty?); "#{tf}[#{s}]" }
      [tf, 'na', tgs] * "\t"
    end
  end

  dep :regulon_modules, :compute => :produce
  dep Viper, :profile, :data => GDSC.gene_expression, :modules => :regulon_modules
  task :gdsc_tf_activity => :tsv do 
    TSV.get_stream(step(:profile))
  end

  dep :regulon_modules
  dep Viper, :profile, :data => CCLE.gene_expression, :modules => :regulon_modules
  task :ccle_tf_activity_Viper => :tsv do
    TSV.get_stream(step(:profile))
  end

  dep :regulon_modules
  dep ROMA, :profile, :data => CCLE.gene_expression, :modules => :regulon_modules
  task :ccle_tf_activity_ROMA => :tsv do
    TSV.get_stream(step(:profile))
  end

  input :tf_inference_method, :select, "Select inference method", "Viper", :select_options => %w(Viper ROMA)
  dep do |jobname, options|
    case options[:tf_inference_method].to_s
    when "Viper"
      {:task => :ccle_tf_activity_Viper, :workflow => CLSS, :inputs => options, :jobname => jobname}
    when "ROMA"
      {:task => :ccle_tf_activity_ROMA, :workflow => CLSS, :inputs => options, :jobname => jobname}
    else
      raise ParameterException, "Unknown method: " << options[:tf_inference_method].inspect
    end
  end
  task :ccle_tf_activity => :tsv do |tf_inference_method|
    TSV.get_stream dependencies.first
  end

  dep :activity
  dep :ccle_tf_activity, :jobname => "Default"
  task :activity_plus_tf => :tsv do |cell_line|
    cell_line = step(:activity).recursive_inputs[:cell_line]
    cl = ccle_cell_line(cell_line)
    raise ParameterException, "Cell line not found: " << cell_line if cl.nil?
    log :cell_line, "Using cell_line #{cl} (#{cell_line})"

    activity = step(:activity).load
    tf_activity = step(:ccle_tf_activity).load

    t_activity = activity.transpose
    t_tf_activity = tf_activity.select(cl).transpose

    t_tf_activity.to_list.each do |gene,v|
      v = v.flatten
      t_activity[gene] = [v]
    end
    t_activity.transpose activity.key_field
  end

  dep :genome
  dep :mRNA
  dep :activity_plus_tf
  dep :st_tf_pathway
  dep Paradigm, :analysis, :genome => :genome, :mRNA => :mRNA, :activity => :activity_plus_tf, :disc => [-0.3, 0.3], :pathway => :st_tf_pathway
  task :steady_states_paradigm => :tsv do
    dumper = TSV::Dumper.new :key_field => "Gene", :fields => %w(Activity), :type => :single
    dumper.init
    TSV.traverse step(:analysis), :into => dumper, :type => :array do |line|
      next if line =~ /^>/
      elem, value = line.split("\t")
      value = value.to_f
      activity = case
                 when value == 0
                   "-"
                 when value < 0
                   "0"
                 when value > 0
                   "1"
                 end

      [elem, activity] 
    end
    dumper
  end

end
