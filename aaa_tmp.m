        cols_abs = length(angles_common2);
        rows = STR_PRE_count+CON_PRE_count+STR_POST_count+CON_POST_count + 1; % adding 1 extra column for holding joint angles
        out_arrays_abs(cols_abs,rows) = zeros;

        for var = 1:length(out_arrays_input_labels)
            
            % reset output arrays
            out_arrays_abs(1:cols_abs,1:rows) = zeros;

            % add as first column, joint angles
            out_arrays_abs(:,1) = angles_common2;
            
            % ROUGH CODING - the following does NOT successfully group
            % subjects into STR/PRE/ETC - subjects are added in the order
            % of the datamaster (same order as all_strength_...)
            
            % add STR PRE subjects first
            for subj = 1:STR_PRE_count
                out_arrays_abs(:,subj+1) = all_strength_isokin_angles(subj,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
            end
            
            % add CON PRE subjects second
            for subj = 1:CON_PRE_count
                out_arrays_abs(:,subj+STR_PRE_count+1) = all_strength_isokin_angles(subj+STR_PRE_count,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
            end

            % add STR POST
            for subj = 1:STR_POST_count
                out_arrays_abs(:,subj+STR_PRE_count+CON_PRE_count+1) = all_strength_isokin_angles(subj+STR_PRE_count+CON_PRE_count,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
            end
            
            % add CON POST
            for subj = 1:CON_POST_count
                out_arrays_abs(:,subj+STR_PRE_count+CON_PRE_count+STR_POST_count+1) = all_strength_isokin_angles(subj+STR_PRE_count+CON_PRE_count+STR_POST_count,out_arrays_input_cols(var):out_arrays_input_cols(var)+14);
            end
            
            % create tables and save as file
            out_arrays_abs_table = array2table(out_arrays_abs,'VariableNames',out_arrays_headers);
            filename_output = strcat('data_output/prism_strength/arrays_strength_acrossangles_', out_arrays_input_labels{var}, '_', datestr(now, 'yyyymmdd_HHMM'));
            writetable(out_arrays_abs_table,filename_output,'Delimiter','\t')

            clear out_arrays_abs_table
        end