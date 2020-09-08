var fileValidated = false;

$(document).one('shiny:connected', function() {
  Shiny.addCustomMessageHandler('querySet', querySet);
  Shiny.addCustomMessageHandler('fileUpload', message => fileValidated = message);
  Shiny.setInputValue('LOAD', true);
  
  // Clickable examples below entries
  loadExamples();
});

$(document).on('shiny:inputchanged', function(event) {
  if(event.name === 'genes')
    validate();
});

// Make the examples clickable
function loadExamples() {
  $('a[genes]').click(function() {
    var sz = $('#genes').selectize()[0].selectize;
    
    sz.clear(true);
    addGenes(JSON.parse($(this).attr('genes')));
  });
}

// Forcibly input genes from a list.
function querySet(genes) { // TODO This is hacky. Might not work
  setTimeout(function() { addGenes(genes); }, 100);
}

// Add genes to the genes bar
function addGenes(genes) {
  var sz = $('#genes').selectize()[0].selectize;
  
  if(genes instanceof Array)
    genes.forEach(gene => sz.createItem(gene, false));
  else
    sz.createItem(genes, false);
      
  Shiny.setInputValue('genes', genes);
}

// Disable the search button if the gene query is empty and no file is uploaded
function validate() {
  var Ng = $('#genes').selectize()[0].selectize.getValue().length;
  var accept = Ng !== 0 || fileValidated;
  
  $('#search').prop('disabled', !accept);
}

function asScientificNotation(row, data) {
  for(var i of [2, 3]) {
    $('td:eq(' + i + ')', row).html(Math.round((data[i] + Number.EPSILON) * 1e3) / 1e3);
  }
}

function onTableCreated() {
  $('[data-toggle="tooltip"]').tooltip();
}

$(document).click(function (e) {
  if($(e.target).parent().find('[data-toggle="popover"]').length > 0) {
    $('[data-toggle="popover"]').popover('hide');
  }
});

function onTableDraw() {
  $('[data-toggle="popover"]').popover();
  
  $('[data-toggle="popover"]').click(function(event) {
    event.preventDefault();
    
    $('[data-toggle="popover"]').not(this).popover('hide');
    $(this).popover('toggle');
  });
  
  $('.spark:not(:has(canvas))').sparkline('html', {
    type: 'bar',
    barColor: '#002145',
    zeroColor: '#002145',
    chartRangeMin: 0,
    chartRangeMax: 1
  });
}

function histogram(data) {
  data = data.split(',').map(x => +x);
  var bins = [];
  
  for(i = 0.1; i <= 1; i += 0.1) {
    bins.push(data.filter(x => x >= (i - 0.1) && x < i).length);
  }
  
  return bins;
}

const mean = arr => arr.reduce((sume, el) => sume + el, 0) / arr.length;

function asSparkline(data, type, row, meta) {
  return type === 'display' ?
    '<span class="spark" style="display: inline-block; width: 100%;" mean="' +
      mean(data.split(',').map(x => +x)) + '">' + histogram(data) + '</span>' :
    data;
}

function unformatSpark(html, row, column, node) {
  if($(node).children().hasClass('spark'))
    return $(html).attr('mean');
  return html;
}

function plotData(e, dt, node, config) {
  Shiny.setInputValue('plotData', {
    page: dt.rows().indexes().toArray(),
    selected: dt.rows('.selected').indexes().toArray()
  }, { priority: 'event' });
}
