<?php

/*
  All Emoncms code is released under the GNU Affero General Public License.
  See COPYRIGHT.txt and LICENSE.txt.

  ---------------------------------------------------------------------
  Emoncms - open source energy visualisation
  Part of the OpenEnergyMonitor project:
  http://openenergymonitor.org
 */

// no direct access
defined('EMONCMS_EXEC') or die('Restricted access');

class PostProcess
{
    private $mysqli;

    public function __construct($mysqli) 
    {
        $this->mysqli = $mysqli;
    }
        
    public function set($userid,$data)
    {
        $userid = (int) $userid;
        // $data = preg_replace('/[^\w\s\-.",:#{}\[\]]/','',$data);
        $data = json_encode($data);
        
        if ($this->get($userid)===false) {
            $stmt = $this->mysqli->prepare("INSERT INTO postprocess ( userid, data ) VALUES (?,?)");
            $stmt->bind_param("is", $userid, $data);
            if ($stmt->execute()) return true; 
        } else {
            $stmt = $this->mysqli->prepare("UPDATE postprocess SET `data`=? WHERE userid=?");
            $stmt->bind_param("si", $data, $userid);
            if ($stmt->execute()) return true;
        }
        return false;
    }
    
    public function get($userid)
    {
        $userid = (int) $userid;
        $result = $this->mysqli->query("SELECT data FROM postprocess WHERE `userid`='$userid'");
        if ($result->num_rows > 0) {
            if ($row = $result->fetch_object()) {
                $data = json_decode($row->data);
                if (!$data || $data==null) $data = array();
                return $data;
            }
        }
        return false;
    }
    
    public function clear_all($userid) 
    {
    
    }

    // Get all processes
    public function get_processes($dir) {

        $processes = array();
    
        $files = scandir($dir);
        for ($i=2; $i<count($files); $i++)
        {
            if (substr($files[$i],-4)==".php" && is_file($dir."/".$files[$i]) && !is_dir($dir."/".$files[$i])) {
                $filename = $files[$i];
                require_once($dir."/".$filename);

                $process_name = str_replace(".php","",$filename);

                if (class_exists("PostProcess_".$process_name)) {
                    $process_class = "PostProcess_".$process_name;
                    $process = new $process_class($dir);

                    if (method_exists($process,"description")) {
                        $process_description = $process->description();
                        if (isset($process_description['settings'])) {
                            $processes[$process_name] = $process_description['settings'];
                        }
                    }
                }
            }
        }
    
        return $processes;
    }
}
