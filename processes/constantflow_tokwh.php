<?php

function constantflow_tokwh($dir,$processitem)
{
    
    $vhc = (float) $processitem->vhc;
    $flow = (float) $processitem->flow;
    $tint = $processitem->tint;
    $text = $processitem->text;
    $out = $processitem->output;
    
    if (!$tint_meta = getmeta($dir,$tint)) return false;
    if (!$text_meta = getmeta($dir,$text)) return false;
    if (!$out_meta = getmeta($dir,$out)) return false;
    
    if (!$tint_fh = @fopen($dir.$tint.".dat", 'rb')) {
        echo "ERROR: could not open $dir $tint.dat\n";
        return false;
    }
    
    if (!$text_fh = @fopen($dir.$text.".dat", 'rb')) {
        echo "ERROR: could not open $dir $text.dat\n";
        return false;
    }
    
    if (!$out_fh = @fopen($dir.$out.".dat", 'c+')) {
        echo "ERROR: could not open $dir $out.dat\n";
        return false;
    }
    
    $compute_meta=compute_meta($tint_meta,$text_meta);
    
    //reading the output meta and if dat file is empty, we adjust interval and start_time
    //we do not report the values in the meta file at this stage. we wait for the dat file to be filled with processed datas
    //if dat file is not empty, meta file should already contain correct values
    print("NOTICE : ouput is : ($out_meta->npoints,$out_meta->interval,$out_meta->start_time) \n");
    if($out_meta->npoints==0) {
        $out_meta->interval=$compute_meta->interval;
        $out_meta->start_time=$compute_meta->start_time;
    }
    print("NOTICE : ouput is now : ($out_meta->npoints,$out_meta->interval,$out_meta->start_time) \n");
    
    $writing_start_time=$out_meta->start_time+($out_meta->interval*$out_meta->npoints);
    $writing_end_time=$compute_meta->writing_end_time;
    $interval=$out_meta->interval;
    $kwh = 0;
    
    if($out_meta->npoints>0) {
        print("NOTICE : file not empty \n");
        $pos_last=($out_meta->npoints-1)*4;
        print("NOTICE : reading last previously processed value in file position $pos_last \n");
        fseek($out_fh,$pos_last);
        $tmp=unpack("f",fread($out_fh,4));
        $kwh = $tmp[1];
        print("NOTICE : accumulation will resume at last previously processed value $kwh \n");
    } else fseek($out_fh,0);
    
    $buffer="";
    
    for ($time=$writing_start_time;$time<$writing_end_time;$time+=$interval){
        
        $pos_tint = floor(($time - $tint_meta->start_time) / $tint_meta->interval);
        $pos_text = floor(($time - $text_meta->start_time) / $text_meta->interval);

        if ($pos_tint>=0 && $pos_tint<$tint_meta->npoints) {
            fseek($tint_fh,$pos_tint*4);
            $tint_tmp = unpack("f",fread($tint_fh,4));
            $value_tint = $tint_tmp[1];
        }
        if ($pos_text>=0 && $pos_text<$text_meta->npoints) {
            fseek($text_fh,$pos_text*4);
            $text_tmp = unpack("f",fread($text_fh,4));
            $value_text = $text_tmp[1];
        }
        
        //if we face one or more NAN value among flow,tint,text, we should not make any calculation and $kwh should remain the same
        if(!is_nan($value_flow) && !is_nan($value_tint) && !is_nan($value_text)){
            $kwh+=0.001*$vhc*$flow*($value_tint-$value_text)*$out_meta->interval/3600;
        }
        //print("$kwh/");
        $buffer.=pack("f",$kwh);
    }
    
    if(!$buffer) {
        print("ERROR: nothing to write - all is up to date \n");
        return false;
    }
    
    if(!$written_bytes=fwrite($out_fh,$buffer)){
        print("ERROR: unable to write to the file with id=$out \n");
        fclose($tint_fh);
        fclose($text_fh);
        fclose($out_fh);
        return false;
    }
    $nbdataswritten=$written_bytes/4;
    print("NOTICE: constantflow_tokwh() wrote $written_bytes bytes ($nbdataswritten float values) \n");
    //we update the meta only as the dat has been filled
    createmeta($dir,$out,$out_meta);
    fclose($tint_fh);
    fclose($text_fh);
    fclose($out_fh);
    print("last time value: $time / $kwh \n");
    updatetimevalue($out,$time,$kwh);
    return true;
}
