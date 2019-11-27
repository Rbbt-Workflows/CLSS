Workflow.require_workflow "Paradigm"
module CLSS

  input :use_tf_activity, :boolean, "Use genomic data", false
  dep :activity
  dep :cl_tf_activity do |jobname,options|
    options[:use_tf_activity] ? {:jobname => jobname, :inputs => options} : []
  end
  task :activity_general => :tsv do |use_tf_activity|
    tsv = step(:activity).load.transpose.to_list
    if use_tf_activity
      tf = step(:cl_tf_activity).load.to_list
      tf.process tf.fields.first do |v|
        case
        when v.to_f  < -1
          -1
        when v.to_f  > 1
          1
        else
          0
        end
      end
                    
      tsv = tsv.annotate(tsv.merge(tf))
    end
    tsv.transpose
  end

  input :use_tf_activity, :boolean, "Use genomic data", false
  input :st_pathway, :text, "Signal transduction pathway in Paradigm"
  dep :st_tf_pathway do |jobname,options|
    if options[:use_tf_activity] 
      {:jobname => jobname, :inputs => options}
    else
      []
    end
  end
  task :paradigm_pathway => :text do |use_tf_activity|
    if use_tf_activity
      step(:st_tf_pathway).load
    else
      pth = self.recursive_inputs[:st_pathway]
      raise RbbtException, "No st_pathway" if pth.nil? or pth.empty?
      pth = pth.read if IO === pth or File === pth
      pt = Misc.is_filename?(pth)? Open.read(pth) : pth
    end
  end

  input :binarize, :boolean, "Binarize output", false
  input :binarize_threshold, :float, "Binarization threshold", 0
  input :binarize_method, :select, "Binarization method", :absolute, :select_options => %w(absolute percentile tails)
  input :use_genome, :boolean, "Use genomic data", false
  input :use_mRNA, :boolean, "Use genomic data", true
  input :use_abundance, :boolean, "Use protein abundance data", true
  input :use_activity, :boolean, "Use activity data", true
  dep :genome do |jobname,options|
    options[:use_genome] ? {:jobname => jobname, :inputs => options} : []
  end
  dep :mRNA do |jobname,options|
    options[:use_mRNA] ? {:jobname => jobname, :inputs => options} : []
  end
  dep :abundance do |jobname,options|
    options[:use_abundance] ? {:jobname => jobname, :inputs => options} : []
  end
  dep :activity_general do |jobname,options|
    options[:use_activity] ? {:jobname => jobname, :inputs => options} : []
  end
  dep :paradigm_pathway, :jobname => 'Default'
  dep Paradigm, :analysis_tsv, :pathway => :paradigm_pathway, :genome => nil, :mRNA => nil, :protein => nil, :activity => nil, :disc => nil, :config_paradigm => nil, :inference => nil do |jobname,options|
    options = options.dup
    IndiferentHash.setup(options)
    options[:genome] = :genome if options[:use_genome]
    options[:mRNA] = :mRNA if options[:use_mRNA]
    options[:protein] = :abundance if options[:use_abundance]
    options[:activity] = :activity_general if options[:use_activity]
    {:inputs => options, :jobname => jobname}
  end
  desc <<-EOF
New version of the steady_states_paradigm

Here we allow the different types of observations to be selected
  EOF
  task :paradigm_ss => :tsv do |binarize,bthreshold,bmethod|
    dumper = TSV::Dumper.new :key_field => "Gene", :fields => %w(Activity), :type => :single, :cast => (binarize ? nil : :to_f)
    dumper.init

    if binarize and bmethod == 'percentile'
      values = []
      TSV.traverse step(:analysis_tsv) do |elem, value|
        value = value.first if Array === value
        value = value.to_f unless value.nil?
        values << value
      end
      l = values.length
      pos = (l * bthreshold / 100).floor
      bthreshold = values.sort[pos]
    end

    if binarize and bmethod == 'tails'
      values = []
      TSV.traverse step(:analysis_tsv) do |elem, value|
        value = value.first if Array === value
        value = value.to_f unless value.nil?
        values << value
      end
      orig = bthreshold
      l = values.length
      pos = (l * orig / 100).floor
      bthreshold = values.sort[pos]

      pos = (l * (1 - orig) / 100).floor
      bthreshold_top = values.sort[pos]
    end

    bthreshold_top ||= bthreshold
    TSV.traverse step(:analysis_tsv), :into => dumper do |elem, value|
      elem = elem.first if Array === elem
      value = value.first if Array === value
      value = value.to_f unless value.nil?
      if binarize
        activity = case
                   when value < bthreshold
                     "0"
                   when value > bthreshold_top
                     "1"
                   else
                     "-"
                   end
      else
        activity = value
      end

      [elem, activity] 
    end
    dumper
  end

  dep :paradigm_ss, :cell_line => :placeholder, :compute => [:bootstrap, nil, :canfail] do |jobname, options|
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
  task :all_paradigm_ss => :tsv do
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
  
  dep :activity, :cell_line => :placeholder, :compute => [:bootstrap, nil, :canfail] do |jobname, options|
    cell_lines = MCLP.RPPA.tsv.fields
    cell_lines.collect do |cell_line|
      {:jobname => cell_line.gsub(' ','_'), :inputs => options.merge(:cell_line => cell_line)}
    end.compact
  end
  task :all_activity => :tsv do
    tsv = nil
    TSV.traverse dependencies, :bar => self.progress_bar("Joining activity files") do |dep|
      next if dep.error?
      cell_line = dep.recursive_inputs[:cell_line]
      this = dep.load.transpose
      this.fields = [cell_line]
      if tsv.nil?
        tsv = this
      else
        tsv.attach this
      end
    end
    good_genes = tsv.keys.select{|k| tsv[k].flatten.reject{|e| e.nil? || e.empty?}.any? }
    tsv.key_field = "Gene"
    tsv.select(good_genes)
  end

  dep :abundance, :cell_line => :placeholder, :compute => [:bootstrap, nil, :canfail] do |jobname, options|
    cell_lines = MCLP.RPPA.tsv.fields
    cell_lines.collect do |cell_line|
      {:jobname => cell_line.gsub(' ','_'), :inputs => options.merge(:cell_line => cell_line)}
    end.compact
  end
  task :all_abundance => :tsv do
    tsv = nil
    TSV.traverse dependencies, :bar => self.progress_bar("Joining activity files") do |dep|
      next if dep.error?
      cell_line = dep.recursive_inputs[:cell_line]
      this = dep.load.transpose
      this.fields = [cell_line]
      if tsv.nil?
        tsv = this
      else
        tsv.attach this
      end
    end
    good_genes = tsv.keys.select{|k| tsv[k].flatten.reject{|e| e.nil? || e.empty?}.any? }
    tsv.key_field = "Gene"
    tsv.select(good_genes)
  end
end
