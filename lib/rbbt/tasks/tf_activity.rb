Workflow.require_workflow "ExTRI"
Workflow.require_workflow "Viper"
Workflow.require_workflow "ROMA"

module CLSS
  dep ExTRI, :regulon, :compute => :produce, :jobname => "Default", :high => true, :skip_post_process => false, :post_process => false, :only_consensus => true, :salvage => false, :tfs => nil
  task :regulon_modules => :text do
    TSV.traverse step(:regulon), :into => :stream do |tf, values|
      tgs = Misc.zip_fields(values).collect{|tf,s|s = 1 if (s.nil? or s.empty?); "#{tf}[#{s}]" }
      [tf, 'na', tgs] * "\t"
    end
  end

  dep :regulon_modules, :compute => :produce
  dep Viper, :profile, :data => RbbtMatrix.new(GDSC.gene_expression).to_name.data_file, :modules => :regulon_modules
  task :gdsc_tf_activity_Viper => :tsv do 
    TSV.get_stream(step(:profile))
  end
  
  dep :regulon_modules
  dep ROMA, :profile, :data => GDSC.gene_expression, :modules => :regulon_modules
  task :gdsc_tf_activity_ROMA => :tsv do
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

  input :cl_data_source, :select, "Source of cell line data", 'CCLE', :select_options => %w(CCLE GDSC)
  input :tf_inference_method, :select, "Select inference method", "Viper", :select_options => %w(Viper ROMA)
  dep :regulon_modules, :compute => :produce
  dep :modules => :regulon_modules do |jobname, options|
    matrix = if options[:cl_data_source] == "GDSC"
               GDSC.gene_expression
             else
               CCLE.gene_expression
             end
    {:workflow => Kernel.const_get(options[:tf_inference_method]), :task => :profile, :jobname => options[:cl_data_source], :inputs => options.merge(:data => matrix)}
  end
  task :global_tf_activity => :tsv do 
    TSV.get_stream(step(:profile))
  end

  input :cell_line, :string, "Cell line name"
  dep :global_tf_activity
  task :cl_tf_activity => :tsv do |cell_line|
    tsv = step(:global_tf_activity).load
    if self.recursive_inputs[:cl_data_source] == "GDSC"
      gdsc_cl = gdsc_cell_line(cell_line)
      raise ParameterException, "GDSC cell line not found: " << cell_line if gdsc_cl.nil?
      log :cell_line, "Using cell_line #{gdsc_cl} (#{cell_line})"
      cosmic_cl = cosmic_cell_line(gdsc_cl)
      raise ParameterException, "COSMIC cell line not found: " << cell_line if cosmic_cl.nil?
      tsv = tsv.transpose.column(cosmic_cl)
      tsv.fields = [cell_line]
    else
      ccle_cl = ccle_cell_line(cell_line)
      raise ParameterException, "CCLE cell line not found: " << cell_line if ccle_cl.nil?
      log :cell_line, "Using cell_line #{ccle_cl} (#{cell_line})"
      tsv = tsv.transpose.column(ccle_cl)
      tsv.fields = [cell_line]
    end
    tsv.cast = :to_f
    tsv
  end

  dep ExTRI, :regulon, :jobname => "Default"
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

  
end
