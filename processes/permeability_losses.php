<?php

chdir("/home/pi/postprocess/lib/");
include "building_model.php";

function permeability_losses($dir,$processitem)
{

    $tint = $processitem->tint;
    $text = $processitem->text;
    $ws = $processitem->ws;
    $out = $processitem->output;
    
    if (!$tint_meta = getmeta($dir,$tint)) return false;
    if (!$text_meta = getmeta($dir,$text)) return false;
    if (!$ws_meta = getmeta($dir,$ws)) return false;
    if (!$out_meta = getmeta($dir,$out)) return false;
    
    if (!$tint_fh = @fopen($dir.$tint.".dat", 'rb')) {
        echo "ERROR: could not open $dir $tint.dat\n";
        return false;
    }
    
    if (!$text_fh = @fopen($dir.$text.".dat", 'rb')) {
        echo "ERROR: could not open $dir $text.dat\n";
        return false;
    }
    
    if (!$ws_fh = @fopen($dir.$ws.".dat", 'rb')) {
        echo "ERROR: could not open $dir $ws.dat\n";
        return false;
    }
    
    if (!$out_fh = @fopen($dir.$out.".dat", 'ab')) {
        echo "ERROR: could not open $dir $out.dat\n";
        return false;
    }
    
    //adjusting the interval of the feed
    $out_interval=max($tint_meta->interval,$text_meta->interval,$ws_meta->interval);
    print("NOTICE : out_interval : $out_interval - original feeds intervals : ($tint_meta->interval,$text_meta->interval,$ws_meta->interval) \n");
    
    //adjusting the start_time of the feed
    $out_start_time=max($tint_meta->start_time,$text_meta->start_time,$ws_meta->start_time);
    $out_start_time=floor($out_start_time/$out_interval)*$out_interval;
    print("NOTICE : out_start_time : $out_start_time - original feeds start_times : ($tint_meta->start_time,$text_meta->start_time,$ws_meta->start_time) \n");
    
    //calculating the supposed end
    $tint_end = $tint_meta->start_time+$tint_meta->npoints*$tint_meta->interval;
    $text_end = $text_meta->start_time+$text_meta->npoints*$text_meta->interval;
    $ws_end = $ws_meta->start_time+$ws_meta->npoints*$ws_meta->interval;
    print("NOTICE : original feeds ends : ($tint_end,$text_end,$ws_end) \n");
    $writing_end_time=min($tint_end,$text_end,$ws_end);
    
    //reading the output meta and if dat file is empty, we adjust interval and start_time
    //we do not report the values in the meta file at this stage. we wait for the dat file to be filled with processed datas
    //if dat file is not empty, meta file should already contain correct values
    print("NOTICE : ouput is : ($out_meta->npoints,$out_meta->interval,$out_meta->start_time) \n");
    if($out_meta->npoints==0) {
        $out_meta->interval=$out_interval;
        $out_meta->start_time=$out_start_time;
    }
    print("NOTICE : ouput is now : ($out_meta->npoints,$out_meta->interval,$out_meta->start_time) \n");
    
    //initializing the model
    $hbat = (float) $processitem->hbat;
    //6 components - towindtop, towindbottom, sidetop, sidebottom, downwindtop, downwindbottom
    //for each component [wind pressure coefficient, component equivalent height in meter, weight of infiltrations on the component, weight of air inlets on the component]
    $components_settings=array(
                [0.25,$hbat/2,1/6,1/6],
                [0.25,0,1/6,1/6],
                [-0.5,$hbat/2,1/6,1/6],
                [-0.5,0,1/6,1/6],
                [-0.5,$hbat/2,1/6,1/6],
                [-0.5,0,1/6,1/6]);
    $q4pasurf = (float) $processitem->q4pasurf;
    $atbat = (float) $processitem->atbat;
    $mea = (float) $processitem->mea;
    
    $home=new Home($q4pasurf,$atbat,$mea,$components_settings);
    
    $qvent = (float) $processitem->qvent;
    
    //calculating the time we are supposed to begin writing
    $writing_start_time=$out_meta->start_time+($out_meta->interval*$out_meta->npoints);
    $interval=$out_meta->interval;
    $buffer="";
    for ($time=$writing_start_time;$time<$writing_end_time;$time+=$interval){
        $pos_tint = floor(($time - $tint_meta->start_time) / $tint_meta->interval);
        $pos_text = floor(($time - $text_meta->start_time) / $text_meta->interval);
        $pos_ws = floor(($time - $ws_meta->start_time) / $ws_meta->interval);
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
        if ($pos_ws>=0 && $pos_ws<$ws_meta->npoints) {
            fseek($ws_fh,$pos_ws*4);
            $ws_tmp = unpack("f",fread($ws_fh,4));
            $value_ws = $ws_tmp[1];
        }
        
        //if we face one or more NAN value among tint,text or ws, we should not make any calculation and just return NAN
        if(!is_nan($value_ws) && !is_nan($value_tint) && !is_nan($value_text)){
            $res=$home->solveur($qvent,$value_ws,$value_tint,$value_text);
            if(is_string($res->zero)) {
                print("ERROR: $res->zero");
                $qvinf=NAN;       
            } else $qvinf=$home->qvinf($res->zero);
        } else $qvinf=NAN;
        //print("$value_ws--$qvinf/");
        $buffer.=pack("f",$qvinf);
    }
    if(!$buffer) {
        print("ERROR: nothing to write - all is up to date \n");
        return false;
    }
    fseek($out_fh,$out_meta->npoints*4);
    if(!$written_bytes=fwrite($out_fh,$buffer)){
        print("ERROR: unable to write to the file with id=$out \n");
        fclose($tint_fh);
        fclose($text_fh);
        fclose($ws_fh);
        fclose($out_fh);
        return false;
    }
    $nbdataswritten=$written_bytes/4;
    print("NOTICE: permeability_losses() wrote $written_bytes bytes ($nbdataswritten float values) \n");
    //we update the meta only as the dat has been filled
    createmeta($dir,$out,$out_meta);
    fclose($tint_fh);
    fclose($text_fh);
    fclose($ws_fh);
    fclose($out_fh);
    print("last time value: $time / $qvinf \n");
    updatetimevalue($out,$time,$qvinf);
    return true;
    
}
 
?>