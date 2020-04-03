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
