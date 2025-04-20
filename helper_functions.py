def get_methylation_percent(df):
    thresh_data = df[["read_id", "call_code"]]
    call_data = thresh_data.groupby(['read_id', 'call_code']).call_code.value_counts().unstack(fill_value=0)
    if 'h' in list(call_data.columns.values):
        call_data['% methylated'] = ((call_data['h'] + call_data['m']) / (call_data['h'] + call_data['m'] + call_data['-'])) * 100
    else:
        call_data['% methylated'] = (call_data['m'] / (call_data['m'] + call_data['-'])) * 100
    call_data = call_data.reset_index()
    return(call_data)

def num_cpgs(df):
    num_cpgs = pd.DataFrame(df['read_id'], columns = ['read_id'])
    num_cpgs['read_length'] = df['read_length']
    num_cpgs = num_cpgs.groupby(num_cpgs.columns.tolist(), as_index = False).size()
    return(num_cpgs)

def read_positions(df):
    start_pos = df['ref_position'] - df['forward_read_position']
    end_pos = (df['read_length'] - df['forward_read_position']) + df['ref_position']
    ref_read_pos = pd.DataFrame(df['read_id'], columns = ['read_id'])
    ref_read_pos["start_pos"] = start_pos
    ref_read_pos["end_pos"] = end_pos
    ref_read_pos["read_length"] = df['read_length']
    ref_read_pos["chrom"] = df['chrom']
    return(ref_read_pos)
    
