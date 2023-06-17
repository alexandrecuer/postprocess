<?php

class PostProcess_trimfeedstart extends PostProcess_common
{
    public function description() {
        return array(
            "name"=>"trimfeedstart",
            "group"=>"Feeds",
            "description"=>"Trim the start of a feed",
            "settings"=>array(
                "feedid"=>array("type"=>"feed", "engine"=>5, "short"=>"Select feed to trim:"),
                "trimtime"=>array("type"=>"value", "short"=>"Enter start time to trim from:")
            )
        );
    }

    public function process($processitem)
    {
        $result = $this->validate($processitem);
        if (!$result["success"]) return $result;
        
        $dir = $this->dir;
        
        $feedid = $processitem->feedid;
        $trimtime = $processitem->trimtime;
        print "TRIM FEED: $feedid, trimtime=$trimtime\n";
        
        $meta = getmeta($dir,$feedid);

        if ($meta->start_time == $trimtime) {
            print "ERROR: feed start_time is already equall to trim time\n";
            return false;
        }
        
        if ($meta->start_time > $trimtime) {
            print "ERROR: feed start_time is more recent than trim time\n";
            return false;
        }

        if (!$if = @fopen($dir.$feedid.".dat", 'rb')) {
            echo "ERROR: could not open $dir$feedid.dat\n";
            return false;
        }
        
        $trimtime = floor($trimtime / $meta->interval) * $meta->interval;
        
        $pos = floor(($trimtime-$meta->start_time) / $meta->interval);
        fseek($if,$pos*4);
        $length = $meta->npoints - $pos;
        
        // Read entire content from trim location forward
        $buffer = fread($if,$length*4.0);
        fclose($if);
        
        if (!$of = @fopen($dir.$feedid.".dat", 'wb')) {
            echo "ERROR: could not open $dir $feedid.dat\n";
            return false;
        }
        
        fwrite($of,$buffer);
        fclose($of);
        
        $meta->start_time = $trimtime;
        $meta->npoints = $length;
        createmeta($dir,$feedid,$meta);
        
        return array("success"=>true);
    }
}
