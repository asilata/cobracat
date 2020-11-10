var data = vis.parseDOTNetwork('digraph ttauta3 {     "abcx" ;     "abcy" ;     "acdy" ;     "acdx" ;     "abdy" ;     "abdx" ;     "bcdy" ;     "bcdx" ;     "abcx" -> "abcx" [taillabel="sa", label="sa:a>b       ", labeltooltip="a:a+x-b;b:b;c:c;x:x     ", url=""];     "abcx" -> "abcx" [taillabel="sa", label="sa:b>a       ", labeltooltip="a:x;b:a;c:c;x:x+b-a     ", url=""];     "abcx" -> "abcx" [taillabel="sb", label="sb:b>c;b-c>x ", labeltooltip="a:a;b:a+b-c-x;c:c;x:x   ", url=""];     "abcx" -> "abcx" [taillabel="sb", label="sb:b>c;x>b-c ", labeltooltip="a:a+x-b+c;b:a;c:c;x:b-c ", url=""];     "abcx" -> "abcx" [taillabel="sc", label="sc:----      ", labeltooltip="b:b,c:b+c+x;a:a;x:x     ", url=""];     "abcx" -> "abcx" [taillabel="ta", label="ta:a>x       ", labeltooltip="b:b;a:b+a-x;c:c;x:x     ", url=""];     "abcx" -> "abcx" [taillabel="ta", label="ta:x>a       ", labeltooltip="b:b+x-a;a:b;c:c;x:a     ", url=""];     "abcx" -> "abcx" [taillabel="tb", label="tb:a>b       ", labeltooltip="c:c;b:c+x;a:b;x:x+a-b   ", url=""];     "abcx" -> "abcx" [taillabel="tb", label="tb:b>a       ", labeltooltip="c:c;b:c+x+b-a;a:a;x:x   ", url=""];     "abcx" -> "abcx" [taillabel="tc", label="tc:c>b;c-b>x ", labeltooltip="c:c-b-x;b:b;a:a;x:x     ", url=""];     "abcx" -> "abcy" [taillabel="sb", label="sb:c>b       ", labeltooltip="a:a+x;b:a;c:b;y:c-b     ", url=""];     "abcx" -> "abdx" [taillabel="tc", label="tc:c>b;x>c-b ", labeltooltip="d:x-c+b;b:b;a:a;x:c-b   ", url=""];     "abcx" -> "abdx" [taillabel="td", label="td:----      ", labeltooltip="a:a;d:a+x;b:b;x:x+c     ", url=""];     "abcx" -> "abdy" [taillabel="tc", label="tc:b>c       ", labeltooltip="d:x;b:c;a:a;y:b-c       ", url=""];     "abcx" -> "bcdy" [taillabel="sd", label="sd:----      ", labeltooltip="c:c+x,d:c;b:b;y:a       ", url=""];     "abcy" -> "abcx" [taillabel="tb", label="tb:a>b       ", labeltooltip="c:c+y;b:c;a:b;x:a-b     ", url=""];     "abcy" -> "abcy" [taillabel="sa", label="sa:a>b;a-b>y ", labeltooltip="a:a-b-y;b:b;c:c;y:y     ", url=""];     "abcy" -> "abcy" [taillabel="sb", label="sb:b>c       ", labeltooltip="a:a;b:a+y+b-c;c:c;y:y   ", url=""];     "abcy" -> "abcy" [taillabel="sb", label="sb:c>b       ", labeltooltip="a:a;b:a+y;c:b;y:y+c-b   ", url=""];     "abcy" -> "abcy" [taillabel="sc", label="sc:c>y       ", labeltooltip="b:b;c:b+c-y;a:a;y:y     ", url=""];     "abcy" -> "abcy" [taillabel="sc", label="sc:y>c       ", labeltooltip="b:b+y-c;c:b;a:a;y:c     ", url=""];     "abcy" -> "abcy" [taillabel="ta", label="ta:----      ", labeltooltip="b:b,a:b+a+y;c:c;y:y     ", url=""];     "abcy" -> "abcy" [taillabel="tb", label="tb:b>a;b-a>y ", labeltooltip="c:c;b:c+b-a-y;a:a;y:y   ", url=""];     "abcy" -> "abcy" [taillabel="tb", label="tb:b>a;y>b-a ", labeltooltip="c:c+y-b+a;b:c;a:a;y:b-a ", url=""];     "abcy" -> "abcy" [taillabel="tc", label="tc:b>c       ", labeltooltip="c:y;b:c;a:a;y:y+b-c     ", url=""];     "abcy" -> "abcy" [taillabel="tc", label="tc:c>b       ", labeltooltip="c:c+y-b;b:b;a:a;y:y     ", url=""];     "abcy" -> "abdx" [taillabel="td", label="td:----      ", labeltooltip="a:a+y,d:a;b:b;x:c       ", url=""];     "abcy" -> "bcdx" [taillabel="sa", label="sa:b>a       ", labeltooltip="d:y;b:a;c:c;x:b-a       ", url=""];     "abcy" -> "bcdy" [taillabel="sa", label="sa:a>b;y>a-b ", labeltooltip="d:y-a+b;b:b;c:c;y:a-b   ", url=""];     "abcy" -> "bcdy" [taillabel="sd", label="sd:----      ", labeltooltip="c:c;d:c+y;b:b;y:y+a     ", url=""];     "abdx" -> "abcx" [taillabel="sc", label="sc:----      ", labeltooltip="b:b;c:b+x;a:a;x:x+d     ", url=""];     "abdx" -> "abcx" [taillabel="sd", label="sd:d>a;x>d-a ", labeltooltip="c:x-d+a;a:a;b:b;x:d-a   ", url=""];     "abdx" -> "abcy" [taillabel="sd", label="sd:a>d       ", labeltooltip="c:x;a:d;b:b;y:a-d       ", url=""];     "abdx" -> "abdx" [taillabel="sa", label="sa:a>b       ", labeltooltip="d:d;a:d+x+a-b;b:b;x:x   ", url=""];     "abdx" -> "abdx" [taillabel="sa", label="sa:b>a       ", labeltooltip="d:d;a:d+x;b:a;x:x+b-a   ", url=""];     "abdx" -> "abdx" [taillabel="sb", label="sb:b>x       ", labeltooltip="a:a;b:a+b-x;d:d;x:x     ", url=""];     "abdx" -> "abdx" [taillabel="sb", label="sb:x>b       ", labeltooltip="a:a+x-b;b:a;d:d;x:b     ", url=""];     "abdx" -> "abdx" [taillabel="sd", label="sd:d>a;d-a>x ", labeltooltip="d:d-a-x;a:a;b:b;x:x     ", url=""];     "abdx" -> "abdx" [taillabel="ta", label="ta:a>d;a-d>x ", labeltooltip="b:b;a:b+a-d-x;d:d;x:x   ", url=""];     "abdx" -> "abdx" [taillabel="ta", label="ta:a>d;x>a-d ", labeltooltip="b:b+x-a+d;a:b;d:d;x:a-d ", url=""];     "abdx" -> "abdx" [taillabel="tb", label="tb:a>b       ", labeltooltip="b:x;a:b;d:d;x:x+a-b     ", url=""];     "abdx" -> "abdx" [taillabel="tb", label="tb:b>a       ", labeltooltip="b:b+x-a;a:a;d:d;x:x     ", url=""];     "abdx" -> "abdx" [taillabel="td", label="td:----      ", labeltooltip="a:a,d:a+d+x;b:b;x:x     ", url=""];     "abdx" -> "abdy" [taillabel="ta", label="ta:d>a       ", labeltooltip="b:b+x;a:b;d:a;y:d-a     ", url=""];     "abdx" -> "acdy" [taillabel="tc", label="tc:----      ", labeltooltip="d:d+x,c:d;a:a;y:b       ", url=""];     "abdy" -> "abcx" [taillabel="sc", label="sc:----      ", labeltooltip="b:b+y,c:b;a:a;x:d       ", url=""];     "abdy" -> "abdx" [taillabel="sa", label="sa:b>a       ", labeltooltip="d:d+y;a:d;b:a;x:b-a     ", url=""];     "abdy" -> "abdy" [taillabel="sa", label="sa:a>b;a-b>y ", labeltooltip="d:d;a:d+a-b-y;b:b;y:y   ", url=""];     "abdy" -> "abdy" [taillabel="sa", label="sa:a>b;y>a-b ", labeltooltip="d:d+y-a+b;a:d;b:b;y:a-b ", url=""];     "abdy" -> "abdy" [taillabel="sb", label="sb:----", labeltooltip="a:a,b:a+b+y;d:d;y:y     ", url=""];     "abdy" -> "abdy" [taillabel="sd", label="sd:a>d       ", labeltooltip="d:y;a:d;b:b;y:y+a-d     ", url=""];     "abdy" -> "abdy" [taillabel="sd", label="sd:d>a       ", labeltooltip="d:d+y-a;a:a;b:b;y:y     ", url=""];     "abdy" -> "abdy" [taillabel="ta", label="ta:a>d       ", labeltooltip="b:b;a:b+y+a-d;d:d;y:y   ", url=""];     "abdy" -> "abdy" [taillabel="ta", label="ta:d>a       ", labeltooltip="b:b;a:b+y;d:a;y:y+d-a   ", url=""];     "abdy" -> "abdy" [taillabel="tb", label="tb:b>a;b-a>y ", labeltooltip="b:b-a-y;a:a;d:d;y:y     ", url=""];     "abdy" -> "abdy" [taillabel="td", label="td:d>y       ", labeltooltip="a:a;d:a+d-y;b:b;y:y     ", url=""];     "abdy" -> "abdy" [taillabel="td", label="td:y>d       ", labeltooltip="a:a+y-d;d:a;b:b;y:d     ", url=""];     "abdy" -> "acdx" [taillabel="tb", label="tb:a>b       ", labeltooltip="c:y;a:b;d:d;x:a-b       ", url=""];     "abdy" -> "acdy" [taillabel="tb", label="tb:b>a;y>b-a ", labeltooltip="c:y-b+a;a:a;d:d;y:b-a   ", url=""];     "abdy" -> "acdy" [taillabel="tc", label="tc:----      ", labeltooltip="d:d;c:d+y;a:a;y:y+b     ", url=""];     "acdx" -> "abdy" [taillabel="sb", label="sb:----", labeltooltip="a:a+x,b:a;d:d;y:c       ", url=""];     "acdx" -> "acdx" [taillabel="sa", label="sa:----      ", labeltooltip="d:d,a:d+a+x;c:c;x:x     ", url=""];     "acdx" -> "acdx" [taillabel="sc", label="sc:c>d       ", labeltooltip="c:c+x-d;d:d;a:a;x:x     ", url=""];     "acdx" -> "acdx" [taillabel="sc", label="sc:d>c       ", labeltooltip="c:x;d:c;a:a;x:x+d-c     ", url=""];     "acdx" -> "acdx" [taillabel="sd", label="sd:d>a;d-a>x ", labeltooltip="c:c;d:c+d-a-x;a:a;x:x   ", url=""];     "acdx" -> "acdx" [taillabel="sd", label="sd:d>a;x>d-a ", labeltooltip="c:c+x-d+a;d:c;a:a;x:d-a ", url=""];     "acdx" -> "acdx" [taillabel="ta", label="ta:a>d;a-d>x ", labeltooltip="a:a-d-x;d:d;c:c;x:x     ", url=""];     "acdx" -> "acdx" [taillabel="tc", label="tc:c>x       ", labeltooltip="d:d;c:d+c-x;a:a;x:x     ", url=""];     "acdx" -> "acdx" [taillabel="tc", label="tc:x>c       ", labeltooltip="d:d+x-c;c:d;a:a;x:c     ", url=""];     "acdx" -> "acdx" [taillabel="td", label="td:c>d       ", labeltooltip="a:a;d:a+x;c:d;x:x+c-d   ", url=""];     "acdx" -> "acdx" [taillabel="td", label="td:d>c       ", labeltooltip="a:a;d:a+x+d-c;c:c;x:x   ", url=""];     "acdx" -> "acdy" [taillabel="sd", label="sd:a>d       ", labeltooltip="c:c+x;d:c;a:d;y:a-d     ", url=""];     "acdx" -> "bcdx" [taillabel="ta", label="ta:a>d;x>a-d ", labeltooltip="b:x-a+d;d:d;c:c;x:a-d   ", url=""];     "acdx" -> "bcdx" [taillabel="tb", label="tb:----      ", labeltooltip="c:c;b:c+x;d:d;x:x+a     ", url=""];     "acdx" -> "bcdy" [taillabel="ta", label="ta:d>a       ", labeltooltip="b:x;d:a;c:c;y:d-a       ", url=""];     "acdy" -> "abdx" [taillabel="sc", label="sc:d>c       ", labeltooltip="b:y;d:c;a:a;x:d-c       ", url=""];     "acdy" -> "abdy" [taillabel="sb", label="sb:----", labeltooltip="a:a;b:a+y;d:d;y:y+c     ", url=""];     "acdy" -> "abdy" [taillabel="sc", label="sc:c>d;y>c-d ", labeltooltip="b:y-c+d;d:d;a:a;y:c-d   ", url=""];     "acdy" -> "acdx" [taillabel="td", label="td:c>d       ", labeltooltip="a:a+y;d:a;c:d;x:c-d     ", url=""];     "acdy" -> "acdy" [taillabel="sa", label="sa:a>y       ", labeltooltip="d:d;a:d+a-y;c:c;y:y     ", url=""];     "acdy" -> "acdy" [taillabel="sa", label="sa:y>a       ", labeltooltip="d:d+y-a;a:d;c:c;y:a     ", url=""];     "acdy" -> "acdy" [taillabel="sc", label="sc:c>d;c-d>y ", labeltooltip="c:c-d-y;d:d;a:a;y:y     ", url=""];     "acdy" -> "acdy" [taillabel="sd", label="sd:a>d       ", labeltooltip="c:c;d:c+y;a:d;y:y+a-d   ", url=""];     "acdy" -> "acdy" [taillabel="sd", label="sd:d>a       ", labeltooltip="c:c;d:c+y+d-a;a:a;y:y   ", url=""];     "acdy" -> "acdy" [taillabel="ta", label="ta:a>d       ", labeltooltip="a:a+y-d;d:d;c:c;y:y     ", url=""];     "acdy" -> "acdy" [taillabel="ta", label="ta:d>a       ", labeltooltip="a:y;d:a;c:c;y:y+d-a     ", url=""];     "acdy" -> "acdy" [taillabel="tc", label="tc:----      ", labeltooltip="d:d,c:d+c+y;a:a;y:y     ", url=""];     "acdy" -> "acdy" [taillabel="td", label="td:d>c;d-c>y ", labeltooltip="a:a;d:a+d-c-y;c:c;y:y   ", url=""];     "acdy" -> "acdy" [taillabel="td", label="td:d>c;y>d-c ", labeltooltip="a:a+y-d+c;d:a;c:c;y:d-c ", url=""];     "acdy" -> "bcdx" [taillabel="tb", label="tb:----      ", labeltooltip="c:c+y,b:c;d:d;x:a       ", url=""];     "bcdx" -> "abcy" [taillabel="ta", label="ta:----      ", labeltooltip="b:b+x,a:b;c:c;y:d       ", url=""];     "bcdx" -> "acdx" [taillabel="sa", label="sa:----      ", labeltooltip="d:d;a:d+x;c:c;x:x+b     ", url=""];     "bcdx" -> "acdx" [taillabel="sb", label="sb:b>c;x>b-c ", labeltooltip="a:x-b+c;c:c;d:d;x:b-c   ", url=""];     "bcdx" -> "acdy" [taillabel="sb", label="sb:c>b       ", labeltooltip="a:x;c:b;d:d;y:c-b       ", url=""];     "bcdx" -> "bcdx" [taillabel="sb", label="sb:b>c;b-c>x ", labeltooltip="b:b-c-x;c:c;d:d;x:x     ", url=""];     "bcdx" -> "bcdx" [taillabel="sc", label="sc:c>d       ", labeltooltip="b:b;c:b+x+c-d;d:d;x:x   ", url=""];     "bcdx" -> "bcdx" [taillabel="sc", label="sc:d>c       ", labeltooltip="b:b;c:b+x;d:c;x:x+d-c   ", url=""];     "bcdx" -> "bcdx" [taillabel="sd", label="sd:d>x       ", labeltooltip="c:c;d:c+d-x;b:b;x:x     ", url=""];     "bcdx" -> "bcdx" [taillabel="sd", label="sd:x>d       ", labeltooltip="c:c+x-d;d:c;b:b;x:d     ", url=""];     "bcdx" -> "bcdx" [taillabel="tb", label="tb:----      ", labeltooltip="c:c,b:c+b+x;d:d;x:x     ", url=""];     "bcdx" -> "bcdx" [taillabel="tc", label="tc:c>b;c-b>x ", labeltooltip="d:d;c:d+c-b-x;b:b;x:x   ", url=""];     "bcdx" -> "bcdx" [taillabel="tc", label="tc:c>b;x>c-b ", labeltooltip="d:d+x-c+b;c:d;b:b;x:c-b ", url=""];     "bcdx" -> "bcdx" [taillabel="td", label="td:c>d       ", labeltooltip="d:x;c:d;b:b;x:x+c-d     ", url=""];     "bcdx" -> "bcdx" [taillabel="td", label="td:d>c       ", labeltooltip="d:d+x-c;c:c;b:b;x:x     ", url=""];     "bcdx" -> "bcdy" [taillabel="tc", label="tc:b>c       ", labeltooltip="d:d+x;c:d;b:c;y:b-c     ", url=""];     "bcdy" -> "abcx" [taillabel="td", label="td:c>d       ", labeltooltip="a:y;c:d;b:b;x:c-d       ", url=""];     "bcdy" -> "abcy" [taillabel="ta", label="ta:----      ", labeltooltip="b:b;a:b+y;c:c;y:y+d     ", url=""];     "bcdy" -> "abcy" [taillabel="td", label="td:d>c;y>d-c ", labeltooltip="a:y-d+c;c:c;b:b;y:d-c   ", url=""];     "bcdy" -> "acdx" [taillabel="sa", label="sa:----      ", labeltooltip="d:d+y,a:d;c:c;x:b       ", url=""];     "bcdy" -> "bcdx" [taillabel="sc", label="sc:d>c       ", labeltooltip="b:b+y;c:b;d:c;x:d-c     ", url=""];     "bcdy" -> "bcdy" [taillabel="sb", label="sb:b>c       ", labeltooltip="b:b+y-c;c:c;d:d;y:y     ", url=""];     "bcdy" -> "bcdy" [taillabel="sb", label="sb:c>b       ", labeltooltip="b:y;c:b;d:d;y:y+c-b     ", url=""];     "bcdy" -> "bcdy" [taillabel="sc", label="sc:c>d;c-d>y ", labeltooltip="b:b;c:b+c-d-y;d:d;y:y   ", url=""];     "bcdy" -> "bcdy" [taillabel="sc", label="sc:c>d;y>c-d ", labeltooltip="b:b+y-c+d;c:b;d:d;y:c-d ", url=""];     "bcdy" -> "bcdy" [taillabel="sd", label="sd:----      ", labeltooltip="c:c,d:c+d+y;b:b;y:y     ", url=""];     "bcdy" -> "bcdy" [taillabel="tb", label="tb:b>y       ", labeltooltip="c:c;b:c+b-y;d:d;y:y     ", url=""];     "bcdy" -> "bcdy" [taillabel="tb", label="tb:y>b       ", labeltooltip="c:c+y-b;b:c;d:d;y:b     ", url=""];     "bcdy" -> "bcdy" [taillabel="tc", label="tc:b>c       ", labeltooltip="d:d;c:d+y;b:c;y:y+b-c   ", url=""];     "bcdy" -> "bcdy" [taillabel="tc", label="tc:c>b       ", labeltooltip="d:d;c:d+y+c-b;b:b;y:y   ", url=""];     "bcdy" -> "bcdy" [taillabel="td", label="td:d>c;d-c>y ", labeltooltip="d:d-c-y;c:c;b:b;y:y     ", url=""]; }');
var allNodes = new vis.DataSet(data.nodes);
var allEdges = new vis.DataSet(data.edges);

var container = document.getElementById("graph");

var nodeoptions = {
    shape: 'box',
    margin: 10,
    font: '14px Helvetica',
};

var edgeoptions = {
    font: {
	size: 0,
	align: 'top',
    },
    smooth: {
	type: 'dynamic'
    },
};

var physicsoptions = {
    enabled: true,
    barnesHut: {
	springLength: 500,
	avoidOverlap: 1,
	springConstant: 0.001
    }
}

var options = {
    autoResize: true,
    height: '100%',
    width: '100%',
    interaction: {
	selectConnectedEdges: false
    },
    nodes: nodeoptions,
    edges: edgeoptions,
    physics: physicsoptions
}

var network;
network = new vis.Network(container, {nodes: allNodes, edges: allEdges}, options);

function undimNode(nodeid){
    return {
	id: nodeid,
	opacity: 1.0,
	font:{
	    color: 'black'
	} 
    }
}

function dimNode(nodeid){
    return {
	id: nodeid,
	opacity: 0.01,
	font:{
	    color: 'rgb(225,225,225)'
	} 
    }
}

function undimEdge(edgeid){
    return {
	id: edgeid,
	font: {
	    size: 14,
	},
	width: 5,
    }
}

function dimEdge(edgeid, width=0){
    return {
	id: edgeid,
	font: {
	    size: 0
	},
	width: width,
    }
}

function dimOrUnDim(params){
    if (params.nodes.length > 0){
	var updateNodeArray = [];
	var keepNodes = Array.from(params.nodes);
	for (n of params.nodes)
	    for (nodeid of network.getConnectedNodes(n))
		keepNodes.push(nodeid);

	//Dim all the nodes except the selected ones
	for (nodeid of allNodes.getIds()){
	    if (keepNodes.includes(nodeid))
		updateNodeArray.push(undimNode(nodeid));
	    else
		updateNodeArray.push(dimNode(nodeid));
	}

	// Dim all the edges except the selected ones and the ones emnating from the selected nodes
	var updateEdgeArray = [];
	var keepEdges = Array.from(params.edges);
	for (nodeid of params.nodes)
	    for (edgeid of network.getConnectedEdges(nodeid))
		keepEdges.push(edgeid)

	for (edgeid of allEdges.getIds()){
	    if (keepEdges.includes(edgeid))
		updateEdgeArray.push(undimEdge(edgeid));
	    else
		updateEdgeArray.push(dimEdge(edgeid));
	}

    }else{
	//Undim all nodes
	var updateNodeArray = [];
	for (nodeid of allNodes.getIds())
	    updateNodeArray.push(undimNode(nodeid));

	//Dim all edges
	var updateEdgeArray = [];
	for (edgeid of allEdges.getIds())
	    updateEdgeArray.push(dimEdge(edgeid, 1));
    }
    allNodes.update(updateNodeArray);
    allEdges.update(updateEdgeArray);    
}

network.on('hold', dimOrUnDim);
