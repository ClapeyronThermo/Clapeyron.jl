abstract type BerthelotRuleModel <: vdW1fRuleModel end

struct BerthelotRule <: EoSParam
end

@newmodelsimple BerthelotRule BerthelotRuleModel BerthelotRule
export vdW1fRule