require 'rbbt/sources/COREAD_phospho_proteome'

module CLSS
  input :cell_line, :string, "Cell line name"
  input :phospho_binarization_type, :select, "Binarization type for phospho levels", :present, :select_options => %w(present threshold_100 levels)
  task :steady_states_roumeliotis => :tsv do |cell_line,type|
    coread_cl = coread_cell_line(cell_line)
    raise ParameterException, "Cell line not found: " << cell_line if coread_cl.nil?
    tsv = TSV.setup({}, "Gene~COREAD Activity#:type=:flat")

    signor = case type.to_sym
             when :present
               COREADPhosphoProteome.signor_activity_present.tsv
             when :threshold_100
               COREADPhosphoProteome.signor_activity_100.tsv
             when :levels
               COREADPhosphoProteome.signor_activity_levels.tsv
             end

    signor.column(coread_cl).each do |site, effect|
      gene = site.split(":").first
      tsv[gene] ||= []
      tsv[gene] << effect
    end
    tsv
  end

  dep :steady_states_roumeliotis
  dep :steady_states_paradigm
  task :compare_coread_paradigm => :tsv do
    sig = step(:steady_states_roumeliotis).load.to_double
    para = step(:steady_states_paradigm).load.to_double

    para = para.attach sig
    para = para.add_field "Match" do |gene, values|
      para_v, coread = values
      
      if coread.empty?
        nil
      elsif para_v.flatten == coread.flatten.uniq
        "MATCH"
      elsif coread.flatten.include? para_v
        "PARTIAL"
      else
        "NO"
      end
    end

    values = para.column("Match").values.flatten.compact
    log :matches, Misc.counts(values).inspect
    set_info :matches, Misc.counts(values)

    para
  end

end
