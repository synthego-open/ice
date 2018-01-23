from ice.classes.sanger_object import SangerObject

def test_trace_to_seq():
    trace = [0,0,1,0]
    trace += [0,0.75, 0.0, 0.25]
    trace += [0.87, 0, 0, 0]
    trace += [0, 0.751, 0.0, 0.249]
    base_order = "ATCG"
    assert "CNAT" == SangerObject.trace_to_base_calls(trace, base_order)