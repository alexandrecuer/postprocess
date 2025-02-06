<?php

/** a model to describe a building - inspired by the THCE CSTB method */

/** sign of a number 
 * @param float $number
 *
 * @return int
 */
function sign( $number ) {
    return ( $number > 0 ) ? 1 : ( ( $number < 0 ) ? -1 : 0 );
}

/** air density in kg/m3 
 * @param float $t air temperature
 * @param float $w air humidity
 *
 * @return float air density in kg/m3
 */
function rhoa($t,$w)
{
    $t=(float)$t;
    $w=0.2;
    $rhoa=(101300-0.37802*$w*exp(59.484085-6790.4985/($t+273.15)-5.02802*log($t+273.15)))/(287.055*($t+273.15));
    return $rhoa;
}
    
/**
 * 
 * home class 
 *
 * model class for a general building - cf THCE method
 *
 * 
*/
class Home
{
    //if you want to use th verbose mode, just bring the two classes Home and Component to a LAMP server. Do not use that mode in EmonCMS
    /** debug variable : 1->verbose, 0->silent */
    CONST VERBOSE=0;
    /** @var float q4pasurf in m3/h/m2 leakage flow rate at differential pressure of 4 pascals divided by atbat - q4pasurf is equivalent to european n50 */
    private $q4pasurf;
    /** @var float atbat wall surface exposed to energy losses in m2*/
    private $atbat;
    /** @var float mea air inlet module in m3/h in m3/h */
    private $mea;
    /** @var array matrix n by 4 of parameters in order to define the n components describing the building. 
    These parameters are : wind pressure coefficient, component equivalent height in meter, weight of infiltrations on the component, weight of air inlets on the component] */
    private $components;
    
    /**
     * constructor
     *
     * @param float $q4pasurf in m3/h/m2 (french n50)
     * @param float $atbat wall surface exposed to energy losses in m2
     * @param float $mea air inlet module in m3/h
     * @param array $components_settings matrix n by 4 
     *
     */
    public function __construct ($q4pasurf,$atbat,$mea,$components_settings)
    {
        //initialize the settings specific to the building
        $this->q4pasurf=(float)$q4pasurf;
        $this->atbat=(float)$atbat;
        $this->mea=(float)$mea;
        
        //initialize all the components from the given array
        for($i=0;$i<sizeof($components_settings);$i++){
            $this->components[$i]=new Component($components_settings[$i]);
        }      
    }
    
    /**
     * calculates permability losses for the whole building ie air fow rate due to infiltrations and air inlets, excluding mechanical ventilation action
     *
     * @param float $pib pressure on the component from inside the building in pa (kg/m/s2)
     * @param float $tint internal temperature in °C
     *
     * @return float permeability losses in m3/h
     */ 
    public function qviea ($pib)
    {
        $q4pasurf=$this->q4pasurf;
        $atbat=$this->atbat;
        $mea=$this->mea;
        $qviea=0;
        for($i=0;$i<sizeof($this->components);$i++) {
            $qvi=$this->components[$i]->_qvi($pib,$q4pasurf,$atbat);
            $qvea=$this->components[$i]->_qvea($pib,$mea);
            //*pow(rhoa($tint,0.2),0.5)*0.8;
            //print("<br>composant".$i.":".$qvi."/".$qvea);
            $qviea+=$qvi+$qvea;
        }
        return $qviea;
    }

    /**
     * calculates infiltration losses for the whole building for the "balance" inside pressure
     *
     * @param float $pib pressure on the component from inside the building in pa (kg/m/s2)
     * @param float $tint internal temperature in °C
     *
     * @return float infiltration losses in m3/h
     */     
    public function qvinf ($pib)
    {
        $q4pasurf=$this->q4pasurf;
        $atbat=$this->atbat;
        $mea=$this->mea;
        $qvinf=0;
        for($i=0;$i<sizeof($this->components);$i++) {
            $qvi=$this->components[$i]->_qvi($pib,$q4pasurf,$atbat);
            $qvea=$this->components[$i]->_qvea($pib,$mea);
            //print("<br>composant".$i.":".$qvi."/".$qvea);
            $qvinf+=max(0,$qvi)+max(0,$qvea);
        }
        return intval($qvinf*10)/10;
    }
    
    /**
     * Calculates internal pressure permitting balance of ventilation and permeability losses, being given tint, text, pext, qvent, known either by calculation or monitoring.
     * Uses a decker-brent algorithm
     *
     * @param float $qvent ventilation losses in m3/h
     * @param float $ws wind speed in m/s
     * @param float $text external temperature in °C
     * @param float $tint internal temperature in °C
     *
     * @return float internal pressure(pa) allowing balance
     */
    public function solveur ($qvent,$ws,$tint,$text)
    {
        $q4pasurf=$this->q4pasurf;
        $atbat=$this->atbat;
        $mea=$this->mea;
        //calculates the external pressure condition on each components and store them
        for($i=0;$i<sizeof($this->components);$i++) {
            $this->components[$i]->_pext($ws,$tint,$text);
            if (self::VERBOSE) $this->components[$i]->show();
        }
        
        //from now we go through a decker-brent engine
        //this algorithm finds the zero of the function f = qvent + qviea
        $a=-1;
        $b=1;
        $res = new stdClass();
        $f_a=$qvent+$this->qviea($a);
        $f_b=$qvent+$this->qviea($b);
        if ($f_a == 0) {
            $res->zero=$a;
            $res->nbit= 0; 
            return $res;
            }
        if ($f_b == 0) {
            $res->zero=$b;
            $res->nbit= 0; 
            return $res;
            }
        //this can happen more than we think with a=-1 and b=1, so we adjust the research interval
        //we suppose the variation direction of the function remains the same
        $offset=25;
        if (sign($f_a)==sign($f_b)) {
            if($f_a<0){
                if ($f_a<$f_b)$b+=$offset;else $a-=$offset;
            }
            if($f_a>0){
                if ($f_a<$f_b)$a-=$offset;else $b+=$offset;
            }
            $f_a=$qvent+$this->qviea($a);
            $f_b=$qvent+$this->qviea($b);
            if (sign($f_a)==sign($f_b)){
                $res->zero="sign of f(a) and f(b) must be opposite tint=$tint, text=$text, ws=$ws \n";
                $res->nbit= 0; 
                return $res;
                }
            }
        if (is_nan($f_a) || is_nan($f_b)) {
            $res->zero="a or b is NaN";
            $res->nbit= 0; 
            return $res;
            }
        $c = $a;
        $f_c = $f_a;
        $d = $b - $c;
        $e = $d;
        $i = 0;
        while ($f_b != 0) {
            if (self::VERBOSE) print("***round $i <br>");
            if (self::VERBOSE) print("next couple abc: ( $a,$b,$c )<br>");
            if (self::VERBOSE) print("f values are: ( $f_a,$f_b,$f_c )<br>"); 
            if (sign($f_a) == sign($f_b)) {
                $a = $c ; $f_a = $f_c;
                $d = $b ; $e = $d;
            }
            if (abs($f_a) < abs($f_b)) {
                $c = $b ; $b = $a ; $a = $c;
                $f_c = $f_b ; $f_b = $f_a ; $f_a = $f_c;
            }
            $m = 0.5*($a - $b);
            if (self::VERBOSE) print("the middle is $m <br>");
            $tol = 2.0*pow(2,-52)*max(abs($b),1.0);
            if ((abs($m) <= $tol) || ($f_b == 0.0)) {
                $res->zero=$b;
                $res->nbit= $i; 
                return $res;
            }

            if (abs($e) < $tol || abs($f_c) <= abs($f_b)) {
                //bisection
                if (self::VERBOSE) print("bisection <br>");
                $d = $m ; 
                $e = $m ;
            } else {
                //interpolation 
                $s = $f_b/$f_c;
                if ($a == $c) {
                    if (self::VERBOSE) print("linear interpolation <br>");
                    $p = 2.0*$m*$s;
                    $q = 1.0 - $s;
                } else {
                    //inverse quadratic interpolation
                    if (self::VERBOSE) print("inverse quadratic interpolation <br>");
                    $q = $f_c / $f_a;
                    $r = $f_b / $f_a;
                    $p = $s*(2.0*$m*$q*($q - $r) - ($b - $c)*($r - 1.0));
                    $q = ($q - 1.0)*($r - 1.0)*($s - 1.0);
                }
                if ($p > 0) $q = -$q; else $p = -$p;
                //is interpolated point acceptable ?
                if ((2.0*$p < 3.0*$m*$q - abs($tol*$q)) && $p < abs(0.5*$e*$q)) {
                    $e = $d;
                    $d = $p/$q;
                } else {
                    $d = $m;
                    $e = $m;
                }
            }
            //next point
            $i = $i + 1;
            $c = $b;
            $f_c = $f_b;
            if (abs($d) > $tol){
                $b = $b + $d;
            } else {
                $b = $b - sign($b-$a)*$tol;
            }
            $f_b = $qvent+$this->qviea($b);
            if ($f_b == 0) {
                $res->zero=$b;
                $res->nbit= $i; 
                return $res;
            }
        }
    }
}

/**
 * 
 * component class 
 *
 * individual brick used to modelize a building
 *
 * 
*/
class Component
{
    /** density of the outside air in kg/m3 */
    CONST RHOAEXT=1.22;
    /** reference temperature of the outside air in K */
    CONST TREF=283;
    /** gravity in m/s2 */
    CONST G=9.81;
    /** @var float wind pressure coefficient */
    private $cp;
    /** @var float component equivalent height in m */
    private $h;
    /** @var float weight of infiltrations on the component */
    private $ri;
    /** @var float weight of air inlets on the component */
    private $rea;
    /** @var float external pressure on the component in pa (kg/m/s2) given windspeed, internal and external temperatures ($ws,$tint,$text)*/
    private $pext;
       
     /**
     * constructor
     *
     * @param array $component_settings [wind pressure coefficient, equivalent height in m, weight of infiltrations, weight of air inlets]
     *
     */
    public function __construct ($component_settings)
    {
        $this->cp=(float)$component_settings[0];
        $this->h=(float)$component_settings[1];
        $this->ri=(float)$component_settings[2];
        $this->rea=(float)$component_settings[3];
    }
    
    /**
     * prints the component static properties except the external pressure
     */     
    public function show(){
        print("<br>");
        print("cp:".$this->cp);
        print(";h:".$this->h);
        print(";ri:".$this->ri);
        print(";rea:".$this->rea);
        print(";pext:".$this->pext);
        print("<br>");       
    }
    
    /**
     * calculates the external pressure on the component
     *
     * @param float $ws wind speed in m/s
     * @param float $tint internal temp in °C or K
     * @param float $text external temp in °C or K
     *
     * @return float the external pressure on the component expressed in pa (kg/m/s2)
     */
    public function _pext($ws,$tint,$text)
    {
        $cp=$this->cp;
        $h=$this->h;
        /** density of the outside air in kg/m3 */
        $rhoaext=self::RHOAEXT;
        /** reference temperature of the outside air in K */
        $tref=self::TREF;
        /** gravity in m/s2 */
        $g=self::G;
        $ws=(intval($ws*10))/10;
        $tint=(intval($tint*10))/10;
        $text=(intval($text*10))/10;
        //$ws=(float)$ws;
        //$tint=(float)$tint;
        //$text=(float)$text;
        $pext = $rhoaext*(0.5*$cp*pow(0.9*$ws,2)-$h*($tint-$text)*$g/$tref);
        //print("/ws:$ws--tint:$tint--text:$text--pext:$pext/");
        $this->pext = $pext;
        return $pext; 
    }
 
    /**
     * calculates the infiltration flow rate
     *
     * @param float $pib pressure on the component from inside the building in pa (kg/m/s2)
     * @param float $q4pasurf q4pasurf in m3/h/m2 - leakage flow rate at differential pressure of 4 pascals divided by atbat - q4pasurf is equivalent to european n50
     * @param float $atbat atbat=wall surface exposed to energy losses in m2
     *
     * @return float infiltration flow rate on the component in m3/h
     */ 
    public function _qvi($pib,$q4pasurf,$atbat)
    {
        $pib=(float)$pib;
        $q4pasurf=(float)$q4pasurf;
        $atbat=(float)$atbat;
        
        $deltap=$this->pext-$pib;
        //print("<br>deltap=".$deltap);
        $qvi=$q4pasurf*$atbat*$this->ri*sign($deltap)*pow(abs($deltap)/4,2/3);
        return $qvi;
       
    }
    
    /**
     * calculates the flow rate of air inlets
     *
     * @param float $pib pressure on the component from inside the building in pa (kg/m/s2)
     * @param float $mea air inlet module in m3/h
     *
     * @return float flow rate of air inlets on the component in m3/h
     */ 
    public function _qvea($pib,$mea)
    {
        $pib=(float)$pib;
        $mea=(float)$mea;
        
        $deltap=$this->pext-$pib;
        //print("/deltap=".$deltap."--pib=".$pib."--pext:".$this->pext."/");
        if ($deltap<20) $qvea=$mea*$this->rea*1.1*sign($deltap)*pow(abs($deltap)/20,0.5);
        if ($deltap>=20) $qvea=$mea*$this->rea*(0.5*$deltap+78)/80;
        return $qvea;
       
    }
    
}

class PostProcess_infiltration_losses extends PostProcess_common
{
    public function description(): array {
        $q4pasurf = "q4pasurf in m3/h/m2
         - leakage flow rate at differential pressure of 4 pascals divided by atbat
         - q4pasurf is equivalent to european n50 :";
        return [
            "name"=>"Infiltration losses",
            "group"=>"Simulation",
            "description"=>"Evaluate permeability losses",
            "settings"=>[
                "tint"=>["type"=>"feed", "engine"=>5, "short"=>"Internal temperature feed :"],
                "text"=>["type"=>"feed", "engine"=>5, "short"=>"External temperature feed :"],
                "ws"=>["type"=>"feed", "engine"=>5, "short"=>"Wind speed feed in m/s :"],
                "qvent"=>["type"=>"value", "short"=>"Ventilation flow rate in m3/h :"],
                "hbat"=>["type"=>"value", "short"=>"building height in m :"],
                "q4pasurf"=>["type"=>"value", "short"=>$q4pasurf],
                "atbat"=>["type"=>"value", "short"=>"atbat in m2 - wall surface exposed to energy losses :"],
                "mea"=>["type"=>"value", "short"=>"mea in m3/h - air inlet module :"],
                "output"=>["type"=>"newfeed", "engine"=>5, "short"=>"Enter output feed name for permeability losses in m3/h :"]
                ]
            ];
        }

    public function infiltration_losses($processitem): array
    {
        $dir = $this->dir;
        $tint = $processitem->tint;
        $text = $processitem->text;
        $ws = $processitem->ws;
        $out = $processitem->output;

        if (!$tint_meta = getmeta($dir,$tint)) return ["success"=>false, "message"=>"could not get meta for $tint"];
        if (!$text_meta = getmeta($dir,$text)) return ["success"=>false, "message"=>"could not get meta for $text"];
        if (!$ws_meta = getmeta($dir,$ws)) return ["success"=>false, "message"=>"could not get meta for $ws"];
        if (!$out_meta = getmeta($dir,$out)) return ["success"=>false, "message"=>"could not get meta for $out"];

        if (!$tint_fh = @fopen($dir.$tint.".dat", 'rb')) {
            return ["success"=>false, "message"=>"could not open $dir $tint.dat"];
        }

        if (!$text_fh = @fopen($dir.$text.".dat", 'rb')) {
            return ["success"=>false, "message"=>"could not open $dir $text.dat"];
        }

        if (!$ws_fh = @fopen($dir.$ws.".dat", 'rb')) {
            return ["success"=>false, "message"=>"could not open $dir $ws.dat"];
        }

        if (!$out_fh = @fopen($dir.$out.".dat", 'ab')) {
            return ["success"=>false, "message"=>"could not open $dir $out.dat"];
        }

        $compute_meta=compute_meta($tint_meta,$text_meta,$ws_meta);

        //reading the output meta and if dat file is empty, we adjust interval and start_time
        //we do not report the values in the meta file at this stage. we wait for the dat file to be filled with processed datas
        //if dat file is not empty, meta file should already contain correct values
        print("NOTICE : ouput is : ($out_meta->npoints,$out_meta->interval,$out_meta->start_time) \n");
        if($out_meta->npoints==0) {
            $out_meta->interval=$compute_meta->interval;
            $out_meta->start_time=$compute_meta->start_time;
        }
        print("NOTICE : ouput is now : ($out_meta->npoints,$out_meta->interval,$out_meta->start_time) \n");

        //initializing the model
        $hbat = (float) $processitem->hbat;
        //6 components - towindtop, towindbottom, sidetop, sidebottom, downwindtop, downwindbottom
        //for each component [wind pressure coefficient, component equivalent height in meter, weight of infiltrations on the component, weight of air inlets on the component]
        $components_settings = [
            [0.25,$hbat/2,1/6,1/6],
            [0.25,0,1/6,1/6],
            [-0.5,$hbat/2,1/6,1/6],
            [-0.5,0,1/6,1/6],
            [-0.5,$hbat/2,1/6,1/6],
            [-0.5,0,1/6,1/6]
        ];
        $q4pasurf = (float) $processitem->q4pasurf;
        $atbat = (float) $processitem->atbat;
        $mea = (float) $processitem->mea;

        $home=new Home($q4pasurf,$atbat,$mea,$components_settings);

        $qvent = (float) $processitem->qvent;

        //calculating the time we are supposed to begin writing
        $writing_start_time=$out_meta->start_time+($out_meta->interval*$out_meta->npoints);
        $writing_end_time=$compute_meta->writing_end_time;
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
            return ["success"=>false, "message"=>"nothing to write - all is up to date"];
        }
        fseek($out_fh,$out_meta->npoints*4);
        if(!$written_bytes=fwrite($out_fh,$buffer)){
            fclose($tint_fh);
            fclose($text_fh);
            fclose($ws_fh);
            fclose($out_fh);
            return ["success"=>false, "message"=>"unable to write to the file with id=$out"];
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
        return [
          "success"=>true,
          "message"=>"bytes written: $written_bytes, last time value: $time, last written value $qvinf"
        ];

    }
}
