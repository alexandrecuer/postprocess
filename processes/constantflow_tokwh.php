<?php

class PostProcess_constantflow_tokwh extends PostProcess_common
{       
    public function description() {
        return array(
            "name"=>"Constant flow to kWh",
            "group"=>"Misc",
            "description"=>"Convert constant flow to kwh",
            "settings"=>array(
                "vhc"=>array("type"=>"value", "short"=>"volumetric heat capacity in Wh/m3/K"),
                "flow"=>array("type"=>"value", "short"=>"constant flow in m3/h"),
                "tint"=>array("type"=>"feed", "engine"=>5, "short"=>"Internal temperature feed / start temperature feed :"),
                "text"=>array("type"=>"feed", "engine"=>5, "short"=>"External temperature feed / return temperature feed :"),
                "output"=>array("type"=>"newfeed", "engine"=>5, "short"=>"Enter output energy feed name (kWh) :")
            )
        );
    }

    public function process($processitem)
    {
        $result = $this->validate($processitem);
        if (!$result["success"]) return $result;
        
        $dir = $this->dir;

        $vhc = (float) $processitem->vhc;
        $flow = (float) $processitem->flow;
        $tint = $processitem->tint;
        $text = $processitem->text;
        $out = $processitem->output;

        if (!$tint_meta = getmeta($dir,$tint)) array("success"=>false, "message"=>"could not open input feed");
        if (!$text_meta = getmeta($dir,$text)) array("success"=>false, "message"=>"could not open input feed");
        if (!$out_meta = getmeta($dir,$out)) array("success"=>false, "message"=>"could not open output feed");

        if (!$tint_fh = @fopen($dir.$tint.".dat", 'rb')) {
            return array("success"=>false, "message"=>"could not open input feed");
        }

        if (!$text_fh = @fopen($dir.$text.".dat", 'rb')) {
            return array("success"=>false, "message"=>"could not open input feed");
        }

        if (!$out_fh = @fopen($dir.$out.".dat", 'c+')) {
            return array("success"=>false, "message"=>"could not open output feed");
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
            $value_tint=NAN;
            $value_text=NAN;
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

            //if we face one or more NAN value among tint,text, we should not make any calculation and $kwh should remain the same
            if(!is_nan($value_tint) && !is_nan($value_text)){
                $kwh+=0.001*$vhc*$flow*max($value_tint-$value_text,0)*$out_meta->interval/3600;
            }
            //print("$kwh/");
            $buffer.=pack("f",$kwh);
        }

        if(!$buffer) {
            return array("success"=>false,"message"=>"nothing to write - all is up to date");
        }

        if(!$written_bytes=fwrite($out_fh,$buffer)){
            fclose($tint_fh);
            fclose($text_fh);
            fclose($out_fh);
            return array("success"=>false,"message"=>"unable to write to the file with id=$out");
        }
        $nbdataswritten=$written_bytes/4;
        print("NOTICE: constantflow_tokwh() wrote $written_bytes bytes ($nbdataswritten float values) \n");
        //we update the meta only as the dat has been filled
        createmeta($dir,$out,$out_meta);
        fclose($tint_fh);
        fclose($text_fh);
        fclose($out_fh);
        print("last time value: $time / $kwh \n");
        
        if ($written_bytes>0) {
            updatetimevalue($out,$time,$kwh);
        }
        return array("success"=>true, "message"=>"bytes written: ".$written_bytes.", last time value: ".$time." ".$kwh);
    }
}