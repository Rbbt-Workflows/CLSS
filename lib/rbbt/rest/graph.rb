class RbbtGraph
  attr_accessor :knowledge_base, :entities, :aesthetics, :associations, :rules

  def initialize
    @entities = {}
    @aesthetics = {}
    @associations = {}
    @rules = {}
  end

  def add_associations(associations, type = :edge)
    @associations[type] ||= []
    @associations[type].concat associations.collect{|i| i }
    @associations[type].uniq
    if AssociationItem === associations
      add_entities associations.target, associations.target_entity_type
      add_entities associations.source, associations.source_entity_type
    end
  end

  def add_entities(entities, type = nil)                                                   
    type = entities.base_entity.to_s if type.nil? and AnnotatedArray === entities
    raise "No type specified and entities are not Annotated, so could not guess" if type.nil? 
    if knowledge_base
      good_entities = knowledge_base.translate(entities, type).compact.uniq
    else
      good_entities = entities
    end
    @namespace ||= entities.organism if entities.respond_to? :organism       
    @entities[type] ||= []
    @entities[type].concat good_entities.to_a
  end                        

  def add_aesthetic(elem, info)
    @aesthetics[elem] ||= []
    @aesthetics[elem] << info
  end

  def add_rule(elem, info)
    @rules[elem] ||= []
    @rules[elem] << info
  end


  def js_model
    js_model = {:entities => {}, :associations => {}, :aes_rules => {}, :edge_aes_rules => {}, :rules => {}, :edge_rules => {}}
    @entities.each do |type, codes|                                                                                          
      js_model[:entities][type] = {:codes => codes, :entity_type => type}
    end
    @associations.each do |type, codes|
      js_model[:associations][type] = {:codes => codes, :database => type}
    end
    @aesthetics.each do |type, info|
      aes_rule_type = (type == :node ? :aes_rules : :edge_aes_rules)
      js_model[aes_rule_type] = info
    end
    @rules.each do |type, info|
      rule_type = (type == :node ? :rules : :edge_rules)
      js_model[rule_type] = info
    end
    js_model
  end
end
