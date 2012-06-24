-- MySQL dump 10.13  Distrib 5.1.63, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: storehouse
-- ------------------------------------------------------
-- Server version	5.1.63-0ubuntu0.11.10.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `ana_in`
--

DROP TABLE IF EXISTS `ana_in`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ana_in` (
  `analysis_id` int(50) DEFAULT NULL,
  `sequence_id` int(50) DEFAULT NULL,
  `unique_id` varchar(100) DEFAULT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `analysis`
--

DROP TABLE IF EXISTS `analysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `analysis` (
  `ID` int(11) NOT NULL AUTO_INCREMENT,
  `input` varchar(50) DEFAULT NULL,
  `db` varchar(50) DEFAULT NULL,
  `program` varchar(50) DEFAULT NULL,
  `parameters` varchar(50) DEFAULT NULL,
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=13 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotation`
--

DROP TABLE IF EXISTS `annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation` (
  `ID` varchar(20) DEFAULT NULL,
  `datatype` varchar(20) DEFAULT NULL,
  `name` varchar(100) DEFAULT NULL,
  `region` varchar(100) DEFAULT NULL,
  `start` int(50) DEFAULT NULL,
  `end` int(50) DEFAULT NULL,
  `right_flank` varchar(20) DEFAULT NULL,
  `rdesc` text,
  `left_flank` varchar(20) DEFAULT NULL,
  `ldesc` text,
  `overlap` varchar(20) DEFAULT NULL,
  `odesc` text,
  `gene_type` varchar(20) DEFAULT NULL,
  `transcript_overlap` varchar(20) DEFAULT NULL,
  `transcript_type` varchar(100) DEFAULT NULL,
  `exon_overlap` varchar(20) DEFAULT NULL,
  `strand` enum('1','-1') DEFAULT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `est`
--

DROP TABLE IF EXISTS `est`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `est` (
  `ID` varchar(10) DEFAULT NULL,
  `name` varchar(100) DEFAULT NULL,
  `datatype` varchar(50) DEFAULT NULL,
  `region` varchar(50) DEFAULT NULL,
  `rstart` int(50) DEFAULT NULL,
  `rend` int(50) DEFAULT NULL,
  `feature` varchar(50) DEFAULT NULL,
  `fstart` int(50) DEFAULT NULL,
  `fend` int(50) DEFAULT NULL,
  `ftype` varchar(50) DEFAULT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_pair`
--

DROP TABLE IF EXISTS `feature_pair`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_pair` (
  `q_id` varchar(100) DEFAULT NULL,
  `h_id` varchar(100) DEFAULT NULL,
  `q_start` int(50) DEFAULT NULL,
  `q_end` int(50) DEFAULT NULL,
  `h_start` int(50) DEFAULT NULL,
  `h_end` int(50) DEFAULT NULL,
  `q_len` int(50) DEFAULT NULL,
  `h_len` int(50) DEFAULT NULL,
  `q_alen` int(50) DEFAULT NULL,
  `h_alen` int(50) DEFAULT NULL,
  `e_value` double DEFAULT NULL,
  `percent_id` float DEFAULT NULL,
  `q_cov` float DEFAULT NULL,
  `h_cov` float DEFAULT NULL,
  `aln` text,
  `strand` enum('1','-1') DEFAULT NULL,
  `frame` enum('-2','-1','0','1','2') DEFAULT NULL,
  `analysis_id` varchar(50) DEFAULT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `rand_feature_pair`
--

DROP TABLE IF EXISTS `rand_feature_pair`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `rand_feature_pair` (
  `q_id` varchar(100) DEFAULT NULL,
  `h_id` varchar(100) DEFAULT NULL,
  `q_start` int(50) DEFAULT NULL,
  `q_end` int(50) DEFAULT NULL,
  `h_start` int(50) DEFAULT NULL,
  `h_end` int(50) DEFAULT NULL,
  `e_value` double DEFAULT NULL,
  `percent_id` int(3) DEFAULT NULL,
  `q_cov` int(3) DEFAULT NULL,
  `h_cov` int(3) DEFAULT NULL,
  `aln` text,
  `strand` enum('1','-1') DEFAULT NULL,
  `frame` enum('-2','-1','0','1','2') DEFAULT NULL,
  `analysis_id` varchar(50) DEFAULT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sequence`
--

DROP TABLE IF EXISTS `sequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sequence` (
  `ID` int(11) NOT NULL AUTO_INCREMENT,
  `dataset` varchar(250) DEFAULT NULL,
  `identifier` varchar(250) DEFAULT NULL,
  `description` varchar(250) DEFAULT NULL,
  `sequence` text,
  `class` varchar(10) DEFAULT NULL,
  PRIMARY KEY (`ID`)
) ENGINE=MyISAM AUTO_INCREMENT=3160665 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2012-06-24 16:51:00
