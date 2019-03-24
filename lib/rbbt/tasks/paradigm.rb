Workflow.require_workflow "Paradigm"
module CLSS


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
      v = [] if v.nil?
      v = v.flatten
      t_activity[gene] = [v]
    end
    t_activity.transpose activity.key_field
  end

  dep :genome
  dep :mRNA
  dep :activity_plus_tf
  dep :st_tf_pathway, :jobname => 'Default'
  dep Paradigm, :analysis, :genome => :genome, :mRNA => :mRNA, :activity => :activity_plus_tf, :disc => [-0.3, 0.3], :pathway => :st_tf_pathway
  input :binarize, :boolean, "Binarize values", true
  task :steady_states_paradigm => :tsv do |binarize|
    dumper = TSV::Dumper.new :key_field => "Gene", :fields => %w(Activity), :type => :single
    dumper.init
    TSV.traverse step(:analysis), :into => dumper, :type => :array do |line|
      next if line =~ /^>/
      elem, value = line.split("\t")
      value = value.to_f
      if binarize
        activity = case
                   when value == 0
                     "-"
                   when value < 0
                     "0"
                   when value > 0
                     "1"
                   end
      else
        activity = value
      end

      [elem, activity] 
    end
    dumper
  end

  dep :steady_states_paradigm, :cell_line => :placeholder, :compute => [:bootstrap, nil, :canfail] do |jobname, options|
    cell_lines = case options[:cl_data_source]
                 when "CCLE"
                   CCLE.cell_lines.keys
                 when "GDSC"
                   GDSC.cell_lines.keys
                 end
    cell_lines.collect do |cell_line|
      {:jobname => cell_line.gsub(' ','_'), :inputs => options.merge(:cell_line => cell_line)}
    end.compact
  end
  task :all_steady_states => :tsv do
    tsv = nil
    dependencies.each do |dep|
      next if dep.error?
      cell_line = dep.recursive_inputs[:cell_line]
      this = dep.load
      this.fields = [cell_line]
      if tsv.nil?
        tsv = this
      else
        tsv.attach this
      end
    end
    tsv
  end

  dep :all_steady_states
  task :all_steady_states_consensus_mayority => :tsv do 
    tsv = step(:all_steady_states).load

    tsv.add_field "Majority vote" do |gene,values|
      num = values.select{|v| v != "-"}.length
      (values.select{|v| v.to_i == 1}.length > num.to_f / 2) ? 1 : 0
    end

    tsv = tsv.slice("Majority vote").to_single
    tsv.fields = ["Activity"]
    tsv
  end
  #{{{ Expresion + genome

  dep :genome
  dep :mRNA
  dep :st_tf_pathway, :jobname => 'Default'
  dep Paradigm, :analysis, :genome => :genome, :mRNA => :mRNA, :activity => nil, :disc => [-0.3, 0.3], :pathway => :st_tf_pathway
  input :binarize, :boolean, "Binarize values", true
  task :steady_states_paradigm_expr => :tsv do |binarize|
    dumper = TSV::Dumper.new :key_field => "Gene", :fields => %w(Activity), :type => :single
    dumper.init
    TSV.traverse step(:analysis), :into => dumper, :type => :array do |line|
      next if line =~ /^>/
      elem, value = line.split("\t")
      value = value.to_f
      if binarize
        activity = case
                   when value == 0
                     "-"
                   when value < 0
                     "0"
                   when value > 0
                     "1"
                   end
      else
        activity = value
      end

      [elem, activity] 
    end
    dumper
  end

  dep :steady_states_paradigm_expr, :cell_line => :placeholder, :compute => [:bootstrap, nil, :canfail] do |jobname, options|
    cell_lines = case options[:cl_data_source]
                 when "CCLE"
                   CCLE.cell_lines.keys
                 when "GDSC"
                   GDSC.cell_lines.keys
                 end
    cell_lines.collect do |cell_line|
      {:jobname => cell_line.gsub(' ','_'), :inputs => options.merge(:cell_line => cell_line)}
    end.compact
  end
  task :all_steady_states_expr => :tsv do
    tsv = nil
    TSV.traverse dependencies, :bar => self.progress_bar("Joining steady state files") do |dep|
      next if dep.error?
      cell_line = dep.recursive_inputs[:cell_line]
      this = dep.load
      this.fields = [cell_line]
      if tsv.nil?
        tsv = this
      else
        tsv.attach this
      end
    end
    tsv
  end

  dep :all_steady_states_expr
  task :all_steady_states_consensus_mayority_expr => :tsv do 
    tsv = step(:all_steady_states_expr).load

    tsv.add_field "Majority vote" do |gene,values|
      num = values.select{|v| v != "-"}.length
      (values.select{|v| v.to_i == 1}.length > num.to_f / 2) ? 1 : 0
    end

    tsv = tsv.slice("Majority vote").to_single
    tsv.fields = ["Activity"]
    tsv
  end


  #{{{ Simple

  dep :mRNA
  dep :st_tf_pathway, :jobname => 'Default'
  dep Paradigm, :analysis, :mRNA => :mRNA, :activity => nil, :disc => [-0.3, 0.3], :pathway => :st_tf_pathway
  input :binarize, :boolean, "Binarize values", true
  task :steady_states_paradigm_simple => :tsv do |binarize|
    dumper = TSV::Dumper.new :key_field => "Gene", :fields => %w(Activity), :type => :single
    dumper.init
    TSV.traverse step(:analysis), :into => dumper, :type => :array do |line|
      next if line =~ /^>/
      elem, value = line.split("\t")
      value = value.to_f
      if binarize
        activity = case
                   when value == 0
                     "-"
                   when value < 0
                     "0"
                   when value > 0
                     "1"
                   end
      else
        activity = value
      end

      [elem, activity] 
    end
    dumper
  end

  dep :steady_states_paradigm_simple, :cell_line => :placeholder, :compute => [:bootstrap, nil, :canfail] do |jobname, options|
    cell_lines = case options[:cl_data_source]
                 when "CCLE"
                   CCLE.cell_lines.keys
                 when "GDSC"
                   GDSC.cell_lines.keys
                 end
    cell_lines.collect do |cell_line|
      {:jobname => cell_line.gsub(' ','_'), :inputs => options.merge(:cell_line => cell_line)}
    end.compact
  end
  task :all_steady_states_simple => :tsv do
    tsv = nil
    TSV.traverse dependencies, :bar => self.progress_bar("Joining steady state files") do |dep|
      next if dep.error?
      cell_line = dep.recursive_inputs[:cell_line]
      this = dep.load
      this.fields = [cell_line]
      if tsv.nil?
        tsv = this
      else
        tsv.attach this
      end
    end
    tsv
  end

  dep :all_steady_states_simple
  task :all_steady_states_consensus_mayority_simple => :tsv do 
    tsv = step(:all_steady_states_simple).load

    tsv.add_field "Majority vote" do |gene,values|
      num = values.select{|v| v != "-"}.length
      (values.select{|v| v.to_i == 1}.length > num.to_f / 2) ? 1 : 0
    end

    tsv = tsv.slice("Majority vote").to_single
    tsv.fields = ["Activity"]
    tsv
  end


  dep :ccle_tf_activity_Viper, :cell_line => :placeholder, :compute => [:bootstrap, 3, :canfail] do |jobname, options|
    CCLE.cell_lines.keys.collect do |cell_line|
      {:jobname => cell_line, :inputs => options.merge(:cell_line => cell_line)}
    end
  end
  task :all_viper_profiles => :tsv do
    tsv = nil
    dependencies.each do |dep|
      next if dep.error?
      cell_line = dep.recursive_inputs[:cell_line]
      this = dep.load
      this.fields = [cell_line]
      if tsv.nil?
        tsv = this
      else
        tsv.attach this
      end
    end
    tsv
  end
  dep :steady_states_paradigm_expr, :cell_line => :placeholder, :compute => [:bootstrap, 3, :canfail] do |jobname, options|
    require 'rbbt/sources/COREAD_phospho_proteome'
    require 'rbbt/ner/rnorm'
    cell_lines = COREADPhosphoProteome.phosphosite_levels.fields
    norm = Normalizer.new(CCLE.cell_lines, :use_keys => true)
    cell_lines.collect do |cell_line|
      begin
        ccle_cl = norm.resolve(cell_line, nil, :max_candidates => 100, :threshold => -100).first
        {:jobname => ccle_cl, :inputs => options.merge(:cell_line => ccle_cl)}
      rescue
        Log.exception $!
      end
    end.compact
  end
  task :CRC_steady_states => :tsv do
    tsv = nil
    dependencies.each do |dep|
      next if dep.error?
      cell_line = dep.recursive_inputs[:cell_line]
      this = dep.load
      this.fields = [cell_line]
      if tsv.nil?
        tsv = this
      else
        tsv.attach this
      end
    end
    tsv
  end

end
