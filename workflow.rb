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

end

require 'rbbt/tasks/data_sources'
require 'rbbt/tasks/tf_activity'
require 'rbbt/tasks/paradigm'
require 'rbbt/tasks/paradigm_new'
require 'rbbt/tasks/expression'
require 'rbbt/tasks/COREAD'

#require 'CLSS/tasks/basic.rb'

#require 'rbbt/knowledge_base/CLSS'
#require 'rbbt/entity/CLSS'

