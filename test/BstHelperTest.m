classdef BstHelperTest < matlab.unittest.TestCase
    
    properties
        tmp_dir
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            tmpd = tempname;
            mkdir(tmpd);
            testCase.tmp_dir = tmpd;
            utest_bst_setup();
        end
    end
    
    methods(TestMethodTeardown)
        function tear_down(testCase)
            rmdir(testCase.tmp_dir, 's');
            utest_clean_bst();
        end
    end
    
    methods(Test)

        function test_run_proc_duplicate_outputs(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            output_name = 'dummy_resampled';
            sFilesOut = bst_process('CallProcess', 'process_resample', sRaw, [], 'freq', 5);
            bst_process('CallProcess', 'process_set_comment', sFilesOut, [], ...
                        'tag', output_name, 'isindex', 0);
            
            sFilesOut = bst_process('CallProcess', 'process_resample', sRaw, [], 'freq', 2);
            bst_process('CallProcess', 'process_set_comment', sFilesOut, [], ...
                        'tag', output_name, 'isindex', 0);
                    
            try        
                [sFilesOut, redone] = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);
            catch ME
                testCase.assertTrue(strcmp(ME.identifier, 'Nst:BstProcOutputError'));
                testCase.assertTrue(~isempty(strfind(ME.message, 'duplicate')));
                testCase.assertTrue(~isempty(strfind(ME.message, 'dummy_raw/dummy_resampled')));
                exception_caught = 1;
            end
            testCase.assertTrue(exception_caught==1);
        end
        
        
        
        function test_run_proc_duplicate_outputs(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            output_name = 'dummy_resampled';
            sFilesOut = bst_process('CallProcess', 'process_resample', sRaw, [], 'freq', 5);
            bst_process('CallProcess', 'process_set_comment', sFilesOut, [], ...
                        'tag', output_name, 'isindex', 0);
            
            sFilesOut = bst_process('CallProcess', 'process_resample', sRaw, [], 'freq', 2);
            bst_process('CallProcess', 'process_set_comment', sFilesOut, [], ...
                        'tag', output_name, 'isindex', 0);
                     
            try
                [sFilesOut, redone] = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);
                throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
            catch ME
                expected_msg = 'Cannot safely manage unique outputs. Found duplicate items: dummy_raw/dummy_resampled';
                testCase.assertTrue(~isempty(strfind(ME.message, expected_msg)));
            end
           
        end
        
        
        
        function test_set_template_default_anat_not_available(testCase)
            %% Ensure that nst_utest protocol exists
            ProtocolName = 'nst_utest';
            % Delete existing protocol
            db_dir = bst_get('BrainstormDbDir');
            gui_brainstorm('DeleteProtocol', ProtocolName);

            nst_protocol_dir = fullfile(db_dir, ProtocolName);
            if exist(nst_protocol_dir, 'dir')
                rmdir(nst_protocol_dir, 's');
            end

            % Create new protocol with default anatomy for all subjects
            gui_brainstorm('CreateProtocol', ProtocolName, 1, 0);
            
            % Make sure template is not available locally
            template_name = 'Dummy_4NIRS';
            template_mri_comment = 'MRI: Dummy 4NIRS';
            template_fn = fullfile(bst_get('BrainstormUserDir'), 'defaults', 'anatomy', [template_name '.zip']);
            if exist(template_fn, 'file')
                delete(template_fn);
            end
            
            nst_bst_set_template_anatomy(template_name, 0, 0);
            
            sSubject = bst_get('Subject', 0);
            testCase.assertMatches(sSubject.Anatomy(1).Comment, template_mri_comment);
        end

        function test_set_template_default_anat_already_available(testCase)
            %% Ensure that nst_utest protocol exists
            ProtocolName = 'nst_utest';
            % Delete existing protocol
            db_dir = bst_get('BrainstormDbDir');
            gui_brainstorm('DeleteProtocol', ProtocolName);

            nst_protocol_dir = fullfile(db_dir, ProtocolName);
            if exist(nst_protocol_dir, 'dir')
                rmdir(nst_protocol_dir, 's');
            end

            % Create new protocol with default anatomy for all subjects
            gui_brainstorm('CreateProtocol', ProtocolName, 1, 0);
            
            % Make sure template is available locally
            template_name = 'Dummy_4NIRS';
            template_bfn = [template_name '.zip'];
            bst_anat_dir = fullfile(bst_get('BrainstormUserDir'), 'defaults', 'anatomy');
            template_fn = fullfile(bst_anat_dir, template_bfn);
            if ~exist(template_fn, 'file')
                template_url = [nst_get_repository_url() '/template/' template_bfn];
                nst_download(template_url, template_fn);
            end
            template_mri_comment = 'MRI: Dummy 4NIRS';            
            nst_bst_set_template_anatomy(template_name, 0, 0);
            
            sSubject = bst_get('Subject', 0);
            testCase.assertMatches(sSubject.Anatomy(1).Comment, template_mri_comment);
        end
 
        function test_set_template_anat(testCase)
            [subject_name, sSubject, iSubject] = bst_create_test_subject('');
            
            template_name = 'Dummy_4NIRS';
            template_mri_comment = 'MRI: Dummy 4NIRS';
            nst_bst_set_template_anatomy(template_name, iSubject, 0);
            
            sSubject = bst_get('Subject', iSubject);
            testCase.assertMatches(sSubject.Anatomy(1).Comment, template_mri_comment);
        end
        
        function test_str_remove_common(testCase)
           
            sl1 = {'common_prefix var1_group common_suffix', ...
                   'common_prefix var2_group common_suffix', ...
                   'common_prefix var3_group common_suffix'};
            
            [strList, commonBegin, commonEnd] = str_remove_common(sl1);
            testCase.assertEqual(length(sl1), length(strList));
            testCase.assertMatches(strList{1}, '1');
            testCase.assertMatches(strList{2}, '2');
            testCase.assertMatches(strList{3}, '3');
            
            testCase.assertMatches(commonBegin, 'common_prefix var');
            testCase.assertMatches(commonEnd, '_group common_suffix');
            
            [strList, commonBegin, commonEnd] = str_remove_common(sl1, 1);
            testCase.assertEqual(length(sl1), length(strList));
            testCase.assertMatches(strList{1}, 'var1');
            testCase.assertMatches(strList{2}, 'var2');
            testCase.assertMatches(strList{3}, 'var3');
            
            testCase.assertMatches(commonBegin, 'common_prefix ');
            testCase.assertMatches(commonEnd, '_group common_suffix');
            
            [strList, commonBegin, commonEnd] = str_remove_common(sl1, 1, ' ');
            testCase.assertEqual(length(sl1), length(strList));
            testCase.assertMatches(strList{1}, 'var1_group');
            testCase.assertMatches(strList{2}, 'var2_group');
            testCase.assertMatches(strList{3}, 'var3_group');
            
            testCase.assertMatches(commonBegin, 'common_prefix ');
            testCase.assertMatches(commonEnd, ' common_suffix');
                
            sl2 = {'var1_group', 'var2_group', 'var3_group'};
            [strList, commonBegin, commonEnd] = str_remove_common(sl2, 1);
            testCase.assertEqual(length(sl2), length(strList));
            testCase.assertMatches(strList{1}, 'var1');
            testCase.assertMatches(strList{2}, 'var2');
            testCase.assertMatches(strList{3}, 'var3');
            
            testCase.assertEmpty(commonBegin);
            testCase.assertMatches(commonEnd, '_group');
            
            [strList, commonBegin, commonEnd] = str_remove_common(sl2, 1, ' ');
            testCase.assertEqual(length(sl2), length(strList));
            testCase.assertMatches(strList{1}, 'var1_group');
            testCase.assertMatches(strList{2}, 'var2_group');
            testCase.assertMatches(strList{3}, 'var3_group');
            
            testCase.assertEmpty(commonBegin, '');
            testCase.assertEmpty(commonEnd, '');  
        end
        
        function test_item_parsing(testCase)
           
           item_str = '';
           item = nst_parse_bst_item_name(item_str);
           
           testCase.assertEmpty(item.condition);
           testCase.assertEmpty(item.subject_name);
           testCase.assertEmpty(item.comment);
 
           item_str = 'data';
           item = nst_parse_bst_item_name(item_str);
           
           testCase.assertEmpty(item.condition);
           testCase.assertEmpty(item.subject_name);
           testCase.assertMatches(item.comment, 'data');

           item_str = 'cond/data';
           item = nst_parse_bst_item_name(item_str);
           
           testCase.assertEmpty(item.subject_name);
           testCase.assertMatches(item.condition, 'cond');
           testCase.assertMatches(item.comment, 'data');
           
           item_str = 'subj/cond/data.mat';
           item = nst_parse_bst_item_name(item_str);
           
           testCase.assertMatches(item.subject_name, 'subj');
           testCase.assertMatches(item.condition, 'cond');
           testCase.assertMatches(item.comment, 'data');

           %TODO: validate error scenarios           
%            item_str = '/waza//data';
%            item = nst_parse_bst_item_name(item_str);
%            
        end
        
        function test_run_proc_force_redo(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            output_name = 'dummy_resampled';
            [sFilesOut, redone] = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);

            testCase.assertTrue(redone==1);
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Set comment');
            
            resampled_data = in_bst_data(sFilesOut);
            testCase.assertMatches(resampled_data.Comment, output_name);
            testCase.assertTrue(all_close(resampled_data.Time, (0:2:(size(y,2)-1))*dt));
            
            %% Call proc again
            [sFilesOut, redone] = nst_run_bst_proc('dummy_resampled', 1, 'process_resample', sRaw, [], 'freq', 5);
            
            testCase.assertTrue(redone==1);
            testCase.assertTrue(exist(file_fullpath(sFilesOut), 'file')==2);
            
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-3,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-3,2}.Comment, 'Delete files');
            testCase.assertTrue(~isempty(strfind(GlobalData.ProcessReports.Reports{end-2,4}, ...
                                                 'Force redo - removed previous result(s):')));
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Set comment');

        end

        
        function test_run_proc_new_cond_subject_dont_redo(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            output_subject = 'dummier';
            output_cond = 'new_cond';
            output_comment = 'dummy_resampled';
            output_name = [output_subject '/' output_cond '/' output_comment];
            [sFilesOut, redone] = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);

            testCase.assertTrue(redone==1);
            
            assert(isempty(nst_get_bst_func_files('test_subject', 'dummy_raw', output_comment)));
            assert(isempty(nst_get_bst_func_files('test_subject', output_cond, output_comment)));
            prev_data_file = nst_get_bst_func_files(output_subject, output_cond, output_comment);
            assert(~isempty(prev_data_file));
            
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-2,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-2,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Set comment');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Move files');
            prev_report_nb_items = length(GlobalData.ProcessReports.Reports);
            
            listing = dir(file_fullpath(prev_data_file));
            prev_data_date = listing.datenum;
            
            resampled_data = in_bst_data(sFilesOut);
            testCase.assertMatches(resampled_data.Comment, output_comment);
            testCase.assertTrue(all_close(resampled_data.Time, (0:2:(size(y,2)-1))*dt));
            
            %% Call proc again
            [sFilesOut, redone] = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);
            
            assert(isempty(nst_get_bst_func_files('test_subject', 'dummy_raw', output_comment)));
            assert(isempty(nst_get_bst_func_files('test_subject', output_cond, output_comment)));
            
            cur_data_file = nst_get_bst_func_files(output_subject, output_cond, output_comment);
            listing = dir(file_fullpath(cur_data_file));
            
            testCase.assertMatches(cur_data_file, prev_data_file);
            testCase.assertEqual(listing.datenum, prev_data_date);
            
            testCase.assertTrue(redone==0);
            
            testCase.assertTrue(exist(file_fullpath(sFilesOut), 'file')==2);
            testCase.assertEqual(length(GlobalData.ProcessReports.Reports), prev_report_nb_items);
        end 
        
        
        function test_head_model_force_redo(testCase)
            %TODO
        end
        
        function test_head_model_redo(testCase)
            global GlobalData;
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]); 
                                    
            %% 1st proc call
            output_name = 'head_model';
            [sFilesOut, redone] = nst_run_bst_proc(output_name, 0, 'process_nst_import_head_model', ...
                                                  sRaw, [], 'use_closest_wl', 1);
                                     
            testCase.assertTrue(redone==1);
            
            testCase.assertEmpty(sFilesOut);
            sInput = bst_process('GetInputStruct', sRaw);
            sStudy = bst_get('Study', sInput.iStudy);
            testCase.assertEqual(sStudy.iHeadModel, 1);
            testCase.assertEqual(sStudy.HeadModel.Comment, output_name);
            
            ilast_report = size(GlobalData.ProcessReports.Reports, 1);
            
             %% Call proc again
             [sFilesOut, redone] = nst_run_bst_proc(output_name, 1, 'process_nst_import_head_model', ...
                                                    sRaw, [], 'use_closest_wl', 1);
             
             testCase.assertTrue(redone==1);
             
             testCase.assertEmpty(sFilesOut);
             testCase.assertTrue(~isempty(strfind(GlobalData.ProcessReports.Reports{ilast_report+1,4}, ...
                                                 'Force redo - removed previous result(s):')));
        end
        
        
        function test_run_proc_dont_redo(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            output_name = 'dummy_resampled';
            [sFilesOut,redone] = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);

            testCase.assertTrue(redone==1);
            
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Set comment');
            prev_report_nb_items = length(GlobalData.ProcessReports.Reports);

            resampled_data = in_bst_data(sFilesOut);
            testCase.assertMatches(resampled_data.Comment, output_name);
            testCase.assertTrue(all_close(resampled_data.Time, (0:2:(size(y,2)-1))*dt));
            
            %% Call proc again
            [sFilesOut, redone] = nst_run_bst_proc('dummy_resampled', 0, 'process_resample', sRaw, [], 'freq', 5);
            
            testCase.assertTrue(redone==0);
            
            testCase.assertTrue(exist(file_fullpath(sFilesOut), 'file')==2);
            testCase.assertEqual(length(GlobalData.ProcessReports.Reports), prev_report_nb_items);
                  
        end        
        
        function test_run_proc_output_mismatch(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            [sFilesOut, redone] = nst_run_bst_proc({'dummy_resampled', 'dummy_dummy'}, 0, 'process_resample', sRaw, [], 'freq', 5);

            testCase.assertTrue(redone==0);
            
            testCase.assertTrue(isempty(sFilesOut));
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Delete files');
            testCase.assertTrue(~isempty(strfind(GlobalData.lastestFullErrMsg, ...
                                'Expected 2 outputs but process produced 1')));                  
        end
    end
end