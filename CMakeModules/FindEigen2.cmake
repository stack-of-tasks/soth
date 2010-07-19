










<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<!-- ViewVC :: http://www.viewvc.org/ -->
<head>
<title>[KDE] Contents of /trunk/KDE/kdeedu/cmake/modules/FindEigen2.cmake</title>
<meta name="generator" content="ViewVC 1.1.6" />
<link rel="shortcut icon" href="/docroot/images/favicon.ico" />
<link rel="stylesheet" href="/docroot/styles.css" type="text/css" />

</head>
<body>
<div class="vc_navheader">
<table><tr>
<td><strong><a href="/?view=roots"><span class="pathdiv">/</span></a><a href="/">[KDE]</a><span class="pathdiv">/</span><a href="/trunk/">trunk</a><span class="pathdiv">/</span><a href="/trunk/KDE/">KDE</a><span class="pathdiv">/</span><a href="/trunk/KDE/kdeedu/">kdeedu</a><span class="pathdiv">/</span><a href="/trunk/KDE/kdeedu/cmake/">cmake</a><span class="pathdiv">/</span><a href="/trunk/KDE/kdeedu/cmake/modules/">modules</a><span class="pathdiv">/</span><a href="/trunk/KDE/kdeedu/cmake/modules/FindEigen2.cmake?view=log">FindEigen2.cmake</a></strong></td>
<td style="text-align: right;"></td>
</tr></table>
</div>
<div style="float: right; padding: 5px;"><a href="http://www.viewvc.org/" title="ViewVC Home"><img src="/docroot/images/viewvc-logo.png" alt="ViewVC logotype" width="240" height="70" /></a></div>
<h1>Contents of /trunk/KDE/kdeedu/cmake/modules/FindEigen2.cmake</h1>

<p style="margin:0;">

<a href="/trunk/KDE/kdeedu/cmake/modules/"><img src="/docroot/images/back_small.png" class="vc_icon" alt="Parent Directory" /> Parent Directory</a>

| <a href="/trunk/KDE/kdeedu/cmake/modules/FindEigen2.cmake?view=log"><img src="/docroot/images/log.png" class="vc_icon" alt="Revision Log" /> Revision Log</a>




</p>

<hr />
<div class="vc_summary">
Revision <a href="/?view=revision&amp;revision=928437"><strong>928437</strong></a> -
(<a href="/trunk/KDE/kdeedu/cmake/modules/FindEigen2.cmake?annotate=928437"><strong>show annotations</strong></a>)
(<a href="/*checkout*/trunk/KDE/kdeedu/cmake/modules/FindEigen2.cmake?revision=928437"><strong>download</strong></a>)


<br /><em>Thu Feb 19 15:43:10 2009 UTC</em>
(16 months, 4 weeks ago)
by <em>ggael</em>







<br />File size: 2571 byte(s)






<pre class="vc_log">check for the correct Eigen2 version
(step will need the next 2.1 version)
</pre>

</div>






<div id="vc_file">
<table cellspacing="0" cellpadding="0">








<tr class="vc_row_odd" id="l1">
<td class="vc_file_line_number">1</td>

<td class="vc_file_line_text"># - Try to find Eigen2 lib
</td>
</tr>




<tr class="vc_row_odd" id="l2">
<td class="vc_file_line_number">2</td>

<td class="vc_file_line_text"># Once done this will define
</td>
</tr>




<tr class="vc_row_odd" id="l3">
<td class="vc_file_line_number">3</td>

<td class="vc_file_line_text">#
</td>
</tr>




<tr class="vc_row_odd" id="l4">
<td class="vc_file_line_number">4</td>

<td class="vc_file_line_text">#  EIGEN2_FOUND - system has eigen lib with correct version
</td>
</tr>




<tr class="vc_row_odd" id="l5">
<td class="vc_file_line_number">5</td>

<td class="vc_file_line_text">#  EIGEN2_INCLUDE_DIR - the eigen include directory
</td>
</tr>




<tr class="vc_row_odd" id="l6">
<td class="vc_file_line_number">6</td>

<td class="vc_file_line_text">#  EIGEN2_VERSION - eigen version
</td>
</tr>




<tr class="vc_row_odd" id="l7">
<td class="vc_file_line_number">7</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l8">
<td class="vc_file_line_number">8</td>

<td class="vc_file_line_text"># Copyright (c) 2006, 2007 Montel Laurent, &lt;montel@kde.org&gt;
</td>
</tr>




<tr class="vc_row_odd" id="l9">
<td class="vc_file_line_number">9</td>

<td class="vc_file_line_text"># Redistribution and use is allowed according to the terms of the BSD license.
</td>
</tr>




<tr class="vc_row_odd" id="l10">
<td class="vc_file_line_number">10</td>

<td class="vc_file_line_text"># For details see the accompanying COPYING-CMAKE-SCRIPTS file.
</td>
</tr>




<tr class="vc_row_odd" id="l11">
<td class="vc_file_line_number">11</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l12">
<td class="vc_file_line_number">12</td>

<td class="vc_file_line_text">include(MacroEnsureVersion)
</td>
</tr>




<tr class="vc_row_odd" id="l13">
<td class="vc_file_line_number">13</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l14">
<td class="vc_file_line_number">14</td>

<td class="vc_file_line_text">if(NOT EIGEN2_MIN_VERSION)
</td>
</tr>




<tr class="vc_row_odd" id="l15">
<td class="vc_file_line_number">15</td>

<td class="vc_file_line_text">  if(NOT Eigen2_FIND_VERSION_MAJOR)
</td>
</tr>




<tr class="vc_row_odd" id="l16">
<td class="vc_file_line_number">16</td>

<td class="vc_file_line_text">    set(Eigen2_FIND_VERSION_MAJOR 2)
</td>
</tr>




<tr class="vc_row_odd" id="l17">
<td class="vc_file_line_number">17</td>

<td class="vc_file_line_text">  endif(NOT Eigen2_FIND_VERSION_MAJOR)
</td>
</tr>




<tr class="vc_row_odd" id="l18">
<td class="vc_file_line_number">18</td>

<td class="vc_file_line_text">  if(NOT Eigen2_FIND_VERSION_MINOR)
</td>
</tr>




<tr class="vc_row_odd" id="l19">
<td class="vc_file_line_number">19</td>

<td class="vc_file_line_text">    set(Eigen2_FIND_VERSION_MINOR 0)
</td>
</tr>




<tr class="vc_row_odd" id="l20">
<td class="vc_file_line_number">20</td>

<td class="vc_file_line_text">  endif(NOT Eigen2_FIND_VERSION_MINOR)
</td>
</tr>




<tr class="vc_row_odd" id="l21">
<td class="vc_file_line_number">21</td>

<td class="vc_file_line_text">  if(NOT Eigen2_FIND_VERSION_PATCH)
</td>
</tr>




<tr class="vc_row_odd" id="l22">
<td class="vc_file_line_number">22</td>

<td class="vc_file_line_text">    set(Eigen2_FIND_VERSION_PATCH 0)
</td>
</tr>




<tr class="vc_row_odd" id="l23">
<td class="vc_file_line_number">23</td>

<td class="vc_file_line_text">  endif(NOT Eigen2_FIND_VERSION_PATCH)
</td>
</tr>




<tr class="vc_row_odd" id="l24">
<td class="vc_file_line_number">24</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l25">
<td class="vc_file_line_number">25</td>

<td class="vc_file_line_text">  set(EIGEN2_MIN_VERSION &quot;${Eigen2_FIND_VERSION_MAJOR}.${Eigen2_FIND_VERSION_MINOR}.${Eigen2_FIND_VERSION_PATCH}&quot;)
</td>
</tr>




<tr class="vc_row_odd" id="l26">
<td class="vc_file_line_number">26</td>

<td class="vc_file_line_text">endif(NOT EIGEN2_MIN_VERSION)
</td>
</tr>




<tr class="vc_row_odd" id="l27">
<td class="vc_file_line_number">27</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l28">
<td class="vc_file_line_number">28</td>

<td class="vc_file_line_text">macro(_eigen2_check_version)
</td>
</tr>




<tr class="vc_row_odd" id="l29">
<td class="vc_file_line_number">29</td>

<td class="vc_file_line_text">  file(READ &quot;${EIGEN2_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h&quot; _eigen2_version_header LIMIT 5000 OFFSET 1000)
</td>
</tr>




<tr class="vc_row_odd" id="l30">
<td class="vc_file_line_number">30</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l31">
<td class="vc_file_line_number">31</td>

<td class="vc_file_line_text">  string(REGEX MATCH &quot;define EIGEN_WORLD_VERSION ([0-9]*)&quot; _eigen2_world_version_match &quot;${_eigen2_version_header}&quot;)
</td>
</tr>




<tr class="vc_row_odd" id="l32">
<td class="vc_file_line_number">32</td>

<td class="vc_file_line_text">  set(EIGEN2_WORLD_VERSION &quot;${CMAKE_MATCH_1}&quot;)
</td>
</tr>




<tr class="vc_row_odd" id="l33">
<td class="vc_file_line_number">33</td>

<td class="vc_file_line_text">  string(REGEX MATCH &quot;define EIGEN_MAJOR_VERSION ([0-9]*)&quot; _eigen2_major_version_match &quot;${_eigen2_version_header}&quot;)
</td>
</tr>




<tr class="vc_row_odd" id="l34">
<td class="vc_file_line_number">34</td>

<td class="vc_file_line_text">  set(EIGEN2_MAJOR_VERSION &quot;${CMAKE_MATCH_1}&quot;)
</td>
</tr>




<tr class="vc_row_odd" id="l35">
<td class="vc_file_line_number">35</td>

<td class="vc_file_line_text">  string(REGEX MATCH &quot;define EIGEN_MINOR_VERSION ([0-9]*)&quot; _eigen2_minor_version_match &quot;${_eigen2_version_header}&quot;)
</td>
</tr>




<tr class="vc_row_odd" id="l36">
<td class="vc_file_line_number">36</td>

<td class="vc_file_line_text">  set(EIGEN2_MINOR_VERSION &quot;${CMAKE_MATCH_1}&quot;)
</td>
</tr>




<tr class="vc_row_odd" id="l37">
<td class="vc_file_line_number">37</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l38">
<td class="vc_file_line_number">38</td>

<td class="vc_file_line_text">  set(EIGEN2_VERSION ${EIGEN2_WORLD_VERSION}.${EIGEN2_MAJOR_VERSION}.${EIGEN2_MINOR_VERSION})
</td>
</tr>




<tr class="vc_row_odd" id="l39">
<td class="vc_file_line_number">39</td>

<td class="vc_file_line_text">  MACRO_ENSURE_VERSION( ${EIGEN2_MIN_VERSION} ${EIGEN2_VERSION} EIGEN2_VERSION_OK)
</td>
</tr>




<tr class="vc_row_odd" id="l40">
<td class="vc_file_line_number">40</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l41">
<td class="vc_file_line_number">41</td>

<td class="vc_file_line_text">  if(NOT EIGEN2_VERSION_OK)
</td>
</tr>




<tr class="vc_row_odd" id="l42">
<td class="vc_file_line_number">42</td>

<td class="vc_file_line_text">  
</td>
</tr>




<tr class="vc_row_odd" id="l43">
<td class="vc_file_line_number">43</td>

<td class="vc_file_line_text">    message(STATUS &quot;Eigen2 version ${EIGEN2_VERSION} found in ${EIGEN2_INCLUDE_DIR}, &quot;
</td>
</tr>




<tr class="vc_row_odd" id="l44">
<td class="vc_file_line_number">44</td>

<td class="vc_file_line_text">                   &quot;but at least version ${EIGEN2_MIN_VERSION} is required&quot;)
</td>
</tr>




<tr class="vc_row_odd" id="l45">
<td class="vc_file_line_number">45</td>

<td class="vc_file_line_text">  endif(NOT EIGEN2_VERSION_OK)
</td>
</tr>




<tr class="vc_row_odd" id="l46">
<td class="vc_file_line_number">46</td>

<td class="vc_file_line_text">endmacro(_eigen2_check_version)
</td>
</tr>




<tr class="vc_row_odd" id="l47">
<td class="vc_file_line_number">47</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l48">
<td class="vc_file_line_number">48</td>

<td class="vc_file_line_text">if (EIGEN2_INCLUDE_DIR)
</td>
</tr>




<tr class="vc_row_odd" id="l49">
<td class="vc_file_line_number">49</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l50">
<td class="vc_file_line_number">50</td>

<td class="vc_file_line_text">  # in cache already
</td>
</tr>




<tr class="vc_row_odd" id="l51">
<td class="vc_file_line_number">51</td>

<td class="vc_file_line_text">  _eigen2_check_version()
</td>
</tr>




<tr class="vc_row_odd" id="l52">
<td class="vc_file_line_number">52</td>

<td class="vc_file_line_text">  set(EIGEN2_FOUND ${EIGEN2_VERSION_OK})
</td>
</tr>




<tr class="vc_row_odd" id="l53">
<td class="vc_file_line_number">53</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l54">
<td class="vc_file_line_number">54</td>

<td class="vc_file_line_text">else (EIGEN2_INCLUDE_DIR)
</td>
</tr>




<tr class="vc_row_odd" id="l55">
<td class="vc_file_line_number">55</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l56">
<td class="vc_file_line_number">56</td>

<td class="vc_file_line_text">find_path(EIGEN2_INCLUDE_DIR NAMES Eigen/Core
</td>
</tr>




<tr class="vc_row_odd" id="l57">
<td class="vc_file_line_number">57</td>

<td class="vc_file_line_text">     PATHS
</td>
</tr>




<tr class="vc_row_odd" id="l58">
<td class="vc_file_line_number">58</td>

<td class="vc_file_line_text">     ${INCLUDE_INSTALL_DIR}
</td>
</tr>




<tr class="vc_row_odd" id="l59">
<td class="vc_file_line_number">59</td>

<td class="vc_file_line_text">     ${KDE4_INCLUDE_DIR}
</td>
</tr>




<tr class="vc_row_odd" id="l60">
<td class="vc_file_line_number">60</td>

<td class="vc_file_line_text">     PATH_SUFFIXES eigen2
</td>
</tr>




<tr class="vc_row_odd" id="l61">
<td class="vc_file_line_number">61</td>

<td class="vc_file_line_text">   )
</td>
</tr>




<tr class="vc_row_odd" id="l62">
<td class="vc_file_line_number">62</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l63">
<td class="vc_file_line_number">63</td>

<td class="vc_file_line_text">if(EIGEN2_INCLUDE_DIR)
</td>
</tr>




<tr class="vc_row_odd" id="l64">
<td class="vc_file_line_number">64</td>

<td class="vc_file_line_text">  _eigen2_check_version()
</td>
</tr>




<tr class="vc_row_odd" id="l65">
<td class="vc_file_line_number">65</td>

<td class="vc_file_line_text">endif(EIGEN2_INCLUDE_DIR)
</td>
</tr>




<tr class="vc_row_odd" id="l66">
<td class="vc_file_line_number">66</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l67">
<td class="vc_file_line_number">67</td>

<td class="vc_file_line_text">include(FindPackageHandleStandardArgs)
</td>
</tr>




<tr class="vc_row_odd" id="l68">
<td class="vc_file_line_number">68</td>

<td class="vc_file_line_text">find_package_handle_standard_args(Eigen2 DEFAULT_MSG EIGEN2_INCLUDE_DIR EIGEN2_VERSION_OK)
</td>
</tr>




<tr class="vc_row_odd" id="l69">
<td class="vc_file_line_number">69</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l70">
<td class="vc_file_line_number">70</td>

<td class="vc_file_line_text">mark_as_advanced(EIGEN2_INCLUDE_DIR)
</td>
</tr>




<tr class="vc_row_odd" id="l71">
<td class="vc_file_line_number">71</td>

<td class="vc_file_line_text">
</td>
</tr>




<tr class="vc_row_odd" id="l72">
<td class="vc_file_line_number">72</td>

<td class="vc_file_line_text">endif(EIGEN2_INCLUDE_DIR)
</td>
</tr>




<tr class="vc_row_odd" id="l73">
<td class="vc_file_line_number">73</td>

<td class="vc_file_line_text">
</td>
</tr>


</table>
</div>





<hr />
<table>
<tr>
<td><address><a href="mailto:<a href="mailto:sysadmin@kde.org">KDE Sysadmin</a>"><a href="mailto:sysadmin@kde.org">KDE Sysadmin</a></a></address></td>
<td style="text-align: right;"><strong><a href="/docroot/help_rootview.html">ViewVC Help</a></strong></td>
</tr>
<tr>
<td>Powered by <a href="http://viewvc.tigris.org/">ViewVC 1.1.6</a></td>
<td style="text-align: right;">&nbsp;</td>
</tr>
</table>
</body>
</html>

